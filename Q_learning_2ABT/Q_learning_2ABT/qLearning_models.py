import torch
import numpy as np

def env_tabt(action_l, state_latent_l: int, outcome_prob_l: float):
    """
    Environment for two-armed bandit task.
    RH/DRH 2022
    
    Args:
        action_l (float):
            Action taken. 0 = right, 1 = left.
        state_latent_l (float):
            Latent state of left arm. 0 = bad, 1 = good.
        outcome_prob_l (float):
            Probability of postive outcome for left arm.

    Returns:
        reward (float):
            Whether positive outcome was received.
        outcome_l (float):
            Whether left arm had positive outcome.
        outcome_r (float):
            Whether right arm had positive outcome.
    """
    outcome_l = torch.rand(1) < ((outcome_prob_l)*(state_latent_l) + (1-outcome_prob_l)*(1-state_latent_l))
    outcome_r = torch.rand(1) < ((1-outcome_prob_l)*(state_latent_l) + (outcome_prob_l)*(1-state_latent_l))
    reward = torch.as_tensor(outcome_l if action_l else outcome_r).type(torch.float32)
    return reward, outcome_l, outcome_r
env_tabt_jit = torch.jit.script(env_tabt)

def sig(x):
    """
    Sigmoid function.
    RH 2022

    Args:
        x (float):
            Input value(s).

    Returns:
        y (float):
            Sigmoid output value(s).
    """
    return 1/(1+torch.exp(-x))
sig_jit = torch.jit.script(sig)

def policy(Q_l, Q_r, action_prev_l, beta, bias_l):
    """
    Q-learning policy for two-armed bandit task.
    RH 2022

    Args:
        Q_l (float):
            Q-value of left arm.
        Q_r (float):
            Q-value of right arm.
        action_prev_l (float):
            Action taken on previous trial. 0 = right, 1 = left.
        beta (float):
            Temperature parameter.
        bias_l (float):
            Bias parameter for left arm.

    Returns:
        prob_l (float):
            'logits': Probability of taking left action.
    """
    return sig_jit(beta*(Q_l - Q_r) + bias_l)
policy_jit = torch.jit.script(policy)

def act(prob_l):
    """
    Action selection. Makes a weighted random choice (boolean)
     based on probability.
    RH 2022

    Args:
        prob_l (float):
            Probability of taking left action.

    Returns:
        action_l (float):
            Action taken. 0 = right, 1 = left.
    """
    action_l = torch.rand(1) < prob_l
    return torch.as_tensor(action_l, dtype=torch.float32)
act_jit = torch.jit.script(act)

def update_Q(reward, action_l, Q_l, Q_r, zeta, alpha, modelID: str='standard'):
    """
    Q-learning update rule.
    RH/DRH 2022

    Args:
        reward (float):
            Whether positive outcome was received.
        action_l (float):
            Action taken. 0 = right, 1 = left.
        Q_l (float):
            Q-value of left arm.
        Q_r (float):
            Q-value of right arm.
        zeta (float):
            discount rate.
        alpha (float):
            learning rate.

    Returns:
        Q_l (float):
            Updated Q-value of left arm.
        Q_r (float):
            Updated Q-value of right arm.
    """

    if modelID == 'standard':
        Q_new_l = (1-(1-action_l)*zeta)*Q_l + alpha*(reward - Q_l)*(action_l)
        Q_new_r = (1-(action_l)*zeta)*Q_r + alpha*(reward - Q_r)*(1-action_l)    
    
    elif modelID == 'reduced':
        Q_new_l = (1-(1-action_l)*alpha)*Q_l + alpha*(reward - Q_l)*(action_l)
        Q_new_r = (1-(action_l)*alpha)*Q_r + alpha*(reward - Q_r)*(1-action_l)
    
    else:
        Q_new_l,Q_new_r = None,None
        raise Exception('not valid model ID')
    

    return Q_new_l, Q_new_r
update_Q_jit = torch.jit.script(update_Q)


def run_session(
    params, 
    mode_generative: bool=False, 
    decisions_emp: torch.Tensor=None, 
    rewards_emp: torch.Tensor=None,
    blockPosition: torch.Tensor=None,
    target: torch.Tensor=None,
    DAB_I_HighProbSel: torch.Tensor=None,
    switch: torch.Tensor=None,
    DAB_I_flipLR_event: torch.Tensor=None,
    state_latent_l=None,
    outcome_prob_l: float=0.8,
    modelID: str='standard',
):
    """
    Run a single session of the two-armed bandit task.
    RH/DRH 2022

    Args:
        params (dict):
            Model parameters.
            Contains:
                beta, bias_l, zeta, alpha
        mode_generative (bool):
            Whether to run in generative mode.
            If True, then state_latent_l and outcome_prob_l
             must be provided.
            If True, decisions_emp and rewards_emp are ignored.
            If True, then will use functions: env_tabt, act
        decisions_emp (float):
            Empirical decisions.
        rewards_emp (float):
            Empirical rewards.
        state_latent_l (array of bool):
            Latent state of left arm. 0 = right, 1 = left.
        outcome_prob_l (float):
            Probability of postive outcome for left arm.
        modelID (str):
            name of model to update Q values

    Returns:
        logger (dict):
            Contains:
                action_l, prob_l, reward, Q_l, Q_r
    """
    n_steps = len(decisions_emp) if not mode_generative else len(state_latent_l)
    logger = {key: torch.ones(n_steps)*torch.nan for key in ['action_l', 'prob_l', 'reward', 'Q_l', 'Q_r', 'blockPosition', 'target', 'DAB_I_HighProbSel', 'switch', 'DAB_I_flipLR_event']}

    Q_l =           torch.tensor(0.)
    Q_r =           torch.tensor(0.)
    action_prev_l = torch.tensor(0.)

    for i_step in range(n_steps):
        ## POLICY, ACTION, REWARD
        prob_l = policy_jit(
            Q_l=Q_l, 
            Q_r=Q_r, 
            action_prev_l=action_prev_l, 
            beta=params[0], 
            bias_l=params[1], 
        )

        if mode_generative:
            action_l = act_jit(prob_l)
            reward, outcome_l, outcome_r = env_tabt_jit(
                action_l=action_l, 
                state_latent_l=state_latent_l[i_step], 
                outcome_prob_l=outcome_prob_l,
            )
        else:
            action_l = decisions_emp[i_step]
            reward = rewards_emp[i_step]

        ## LOG TRIAL VARIABLES
        logger['action_l'][i_step] = action_l
        logger['prob_l'][i_step] = prob_l
        logger['reward'][i_step] = reward
        logger['Q_l'][i_step] = Q_l
        logger['Q_r'][i_step] = Q_r
        logger['blockPosition'][i_step] = blockPosition[i_step]
        logger['target'][i_step] = target[i_step]
        logger['DAB_I_HighProbSel'][i_step] = DAB_I_HighProbSel[i_step]
        logger['switch'][i_step] = switch[i_step]
        logger['DAB_I_flipLR_event'][i_step] = DAB_I_flipLR_event[i_step]
        ## UPDATE Q VALUES
        Q_l, Q_r = update_Q_jit(
            reward=reward, 
            action_l=action_l, 
            Q_l=Q_l, 
            Q_r=Q_r, 
            zeta=params[2], 
            alpha=params[3],
            modelID=modelID
        )
        action_prev_l = action_l
    return logger
run_session_jit = torch.jit.script(run_session)


def epoch_step(optimizer, fn_loss, loss_rolling, params, decisions_emp, rewards_emp, modelID):
    """
    Run a single epoch of the model. This calls run_session once and optimizes
     the model parameters.
    RH 2022

    Args:
        optimizer (torch.optim):
            Optimizer object.
            Should be similar to torch.optim.SGD
        fn_loss (function):
            Loss function.
            Should be similar to torch.nn.CrossEntropyLoss
        loss_rolling (float):
            List of rolling loss values.
        params (dict):
            Model parameters.
            See run_session for details on this arg.
        decisions_emp (float):
            Empirical decisions.
        rewards_emp (float):
            Empirical rewards.
        temp_CEL (float):
            Temperature parameter for cross-entropy loss.

    Returns:
        loss_rolling (float):
            List of rolling loss values.
    """
    optimizer.zero_grad()
    # logger = run_session(decisions_emp, rewards_emp, params)
    logger = run_session_jit(
        params=params, 
        mode_generative=False, 
        decisions_emp=decisions_emp, 
        rewards_emp=rewards_emp,
        state_latent_l=None,
        outcome_prob_l=0.8,
        modelID = modelID
    )

    probs = torch.vstack([1-logger['prob_l'], logger['prob_l']]).T
    # loss = fn_loss(probs/temp_CEL, decisions_emp.type(torch.long))
    loss = fn_loss(probs.log(), decisions_emp.type(torch.long))
    loss.backward()
    optimizer.step()
    
    loss_rolling.append(loss.item())
    return logger, loss_rolling
# epoch_step_jit = torch.jit.script(epoch_step)

def epoch_step_batch(
        optimizer, 
        fn_loss, 
        loss_rolling, 
        params, 
        decisions_emp, 
        rewards_emp, 
        blockPosition,
        target,
        DAB_I_HighProbSel,
        switch,
        DAB_I_flipLR_event,
        modelID='standard',
    ):
    """
    Run a single epoch of the model. This calls run_session once and optimizes
     the model parameters.
    RH 2022

    Args:
        optimizer (torch.optim):
            Optimizer object.
            Should be similar to torch.optim.SGD
        fn_loss (function):
            Loss function.
            Should be similar to torch.nn.CrossEntropyLoss
        loss_rolling (float):
            List of rolling loss values.
        params (dict):
            Model parameters.
            See run_session for details on this arg.
        decisions_emp (list of float):
            List of empirical decision arrays.
        rewards_emp (list of float):
            List of empirical reward arrays.
        blockPosition (dataframe or dict of arrays):
            Dataframe of values of same length as decisions_emp
             and rewards_emp  (len is shape[0])
        modelID (str):
            Type of model to run.

    Returns:
        loss_rolling (float):
            List of rolling loss values.
    """
    assert len(decisions_emp) == len(rewards_emp), "RH ERROR: decisions_emp and rewards_emp must be same length"
    assert len(blockPosition) == len(decisions_emp), "RH ERROR: blockPosition should be an array with first dimension length equal to len(decisions_emp)"

    optimizer.zero_grad()
    logger = [run_session_jit(params=params, decisions_emp=d, rewards_emp=r, blockPosition=b, target=t, DAB_I_HighProbSel=h, switch =s, DAB_I_flipLR_event=f, modelID=modelID) for d, r, b, t, h, s, f in zip(decisions_emp, rewards_emp, blockPosition, target, DAB_I_HighProbSel, switch, DAB_I_flipLR_event)]
    probs = [torch.vstack([1-l['prob_l'], l['prob_l']]).T for l in logger]
    loss = fn_loss(torch.vstack(probs).log(), torch.cat(decisions_emp).type(torch.long))
    loss.backward()
    optimizer.step()

    loss_rolling.append(loss.item())
    return logger, loss_rolling

class Convergence_checker:
    """
    'Convergence Checker' Class.
    Checks for convergence of the optimization. Uses Ordinary Least Squares (OLS) to 
     fit a line to the last 'window_convergence' number of iterations.
    """
    def __init__(
        self,
        tol_convergence=1e-2,
        window_convergence=100,
    ):
        """
        Initialize the convergence checker.

        Args:
            tol_convergence (float): 
                Tolerance for convergence.
                Corresponds to the slope of the line that is fit.
            window_convergence (int):
                Number of iterations to use for fitting the line.
        """
        self.window_convergence = window_convergence
        self.tol_convergence = tol_convergence

        self.line_regressor = torch.cat((torch.linspace(0,1,window_convergence)[:,None], torch.ones((window_convergence,1))), dim=1)

    def OLS(self, y):
        """
        Ordinary least squares.
        Fits a line to y.
        """
        X = self.line_regressor
        theta = torch.inverse(X.T @ X) @ X.T @ y
        y_rec = X @ theta
        bias = theta[-1]
        theta = theta[:-1]

        return theta, y_rec, bias

    def __call__(
        self,
        loss_history,
    ):
        """
        Forward pass of the convergence checker.
        Checks if the last 'window_convergence' number of iterations are
         within 'tol_convergence' of the line fit.

        Args:
            loss_history (list):
                List of loss values for the last 'window_convergence' number of iterations.

        Returns:
            diff_window_convergence (float):
                Difference of the fit line over the range of 'window_convergence'.
            loss_smooth (float):
                The mean loss over 'window_convergence'.
            converged (bool):
                True if the 'diff_window_convergence' is less than 'tol_convergence'.
        """
        if len(loss_history) < self.window_convergence:
            return np.nan, np.nan, False
        loss_window = torch.as_tensor(loss_history[-self.window_convergence:], device='cpu', dtype=torch.float32)
        theta, y_rec, bias = self.OLS(y=loss_window)

        diff_window_convergence = (y_rec[-1] - y_rec[0])
        loss_smooth = loss_window.mean()
        converged = True if torch.abs(diff_window_convergence) < self.tol_convergence else False
        return diff_window_convergence.item(), loss_smooth.item(), converged

