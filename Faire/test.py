import os,sys
from hmmpytk import hmm_faster

hmm_model = hmm_faster.HMM()
hmm_model.set_states(['NN', 'NP', 'VBD', 'JJ'])
hmm_model.set_observations(['lamp', 'light', 'race', 'huge'])
hmm_model.randomize_matrices(seed = 100)
Pi_matrix={'NN': 1, 'NP': 0, 'VBD': 0.1, 'JJ':0.3}
hmm_model.set_initial_matrix(Pi_matrix)
print hmm_model.viterbi(['lamp', 'light', 'race', 'huge'])
print hmm_model.get_model()
hmm_model.train(['light', 'race', 'huge','light','race'], max_iteration=100, delta=0.0001)
print hmm_model.get_model()
#hmm_model.viterbi(['them', ',', 'they', 'asked', 'him', 'no', 'question', ',', 'for', 'his', 'face', 'told', 'them', 'everything', '.'])
#hmm_model.train(['PRP', ',', 'PRP', 'VBD', 'PRP', 'DT', 'NN', ',', 'IN', 'PRP$', 'NN', 'VBD', 'PRP', 'NN', '.' ], max_iteration=100, delta=0.0001)
