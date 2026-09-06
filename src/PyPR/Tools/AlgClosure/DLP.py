from PyPR.FeedbackFunctions import MPR
from PyPR import FeedbackRegister

import numpy as np

def DLP_brute_force(
    primitive_polynomial,
    base_element,
    limit = 2**30
):
    def dlp(search_element):
        search_element = np.array(search_element,dtype='uint8')
        M = MPR(
            len(primitive_polynomial)-1,
            primitive_polynomial,
            base_element
        )

        M.compile()
        F = FeedbackRegister(1,M)
        for i,state in enumerate(F.run(limit)):
            if np.all(state._state == search_element):
                return i
        return None
    return dlp

print(DLP_brute_force(
    [1,1,0,1],
    [0,1,0]
)([1,1,1]))