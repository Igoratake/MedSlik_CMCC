import numpy as np
import os
from BayesianOptimizationLib_.bayes_opt import BayesianOptimization
from util_lib_ import find_max
from BayesianOptimizationLib_.bayes_opt.path_ import Path
from multi_fss_processing_cor_ import compute_fss

def workflow(parent):
    
    ''' 
    Function that executes the launch of the model with 
    consequent extraction of the fss related to the run 
    '''
    def medslikII(wind_drag, wind_angle, horizontal_diff):

        params = {
            'wind_drag': np.float32(wind_drag),
            'wind_angle': np.float32(wind_angle),
            'horizontal_diff': np.float32(horizontal_diff)
        }

        ''' 
        MEDSLIK II
        Model launcher for try new parameter returned from Bayesian Optimization
        '''
        os.system(Path.LAUNCH_MEDSLIKII)

        return compute_fss(parent)

    '''
    Definition of the initial bounds for performing the bayesian optimization
    '''
    bounds = {
        # 'delt':(17,94), non modificabile
        'wind_drag' : (0, 0.05),
        'wind_angle' : (0, 15), #(-15,15), for northern hemisphere
        'horizontal_diff' : (0, 20)
    }

    '''
    Definition of the optimizer that performs the Bayesian Optimization
    '''
    optimizer = BayesianOptimization(
        f=medslikII,
        pbounds=bounds,
        random_state=10,
        verbose=2
    )

    '''
    Launch the maximize function who try to find better parameters for the simulation
    '''
    optimizer.maximize(init_points=1, n_iter=1)

    '''
    Return the best result found by Bayesian Optimization
    '''
    find_max(Path.WORKFLOW_FOLDER + Path.FINAL_FSS_RESULT_FILE)

    #optimizer.max