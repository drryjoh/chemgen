
# import math
# import numpy as np
# import matplotlib.pyplot as plt
# import scipy as sp
# import tensorflow as tf
# from tensorflow import keras
# # from mlp import custom_loss  


# # load model and input matrix
# model = keras.models.load_model('../mlp_inv_model_norm.keras')
# model.summary()
# A = np.loadtxt('../data_200/A_0.csv', delimiter=',').astype(np.float32)
# A = A.reshape(-1, 96, 96, 1)

# # perform prediction
# predictions = model.predict(A)
# # np.set_printoptions(threshold=np.inf)
# print(predictions)



























































# #!/usr/bin/env python3
# import numpy as np
# from tensorflow import keras
# import os

# MODEL_PATH     = '../MLP_3.keras'
# STATS_CSV_PATH = '../MLP_3.csv'
# A_PATH         = '../data_200/A_199.csv'
# EPS = 1e-8
# OUTPUT_FOLDER = "./py_layer_outputs"

# # input model
# model = keras.models.load_model(MODEL_PATH)

# # load csv filess
# X_mean = X_std = y_mean = y_std = None
# with open(STATS_CSV_PATH, 'r') as f:
#     for line in f:
#         # seperate key and value
#         key, rest = line.split(':', 1)
#         # seperate delimiters
#         vals_str = rest.strip().lstrip('[').rstrip(']')
#         arr = np.fromstring(vals_str, sep=',', dtype=np.float32)
#         if key.strip() == 'input_mean':
#             X_mean = arr
#         elif key.strip() == 'input_std':
#             X_std = arr
#         elif key.strip() == 'output_mean':
#             y_mean = arr
#         elif key.strip() == 'output_std':
#             y_std = arr

# # sanity check
# for name, v in [('input_mean',X_mean),('input_std',X_std),('output_mean',y_mean),('output_std',y_std)]:
#     assert v is not None, f"{name} not found in {STATS_CSV_PATH}"

# # load data
# A = np.loadtxt(A_PATH, delimiter=',', dtype=np.float32) 
# A_flat = A.ravel()                     

# # normalize data
# A_norm = (A_flat - X_mean) / (X_std + EPS)
# A_norm = np.expand_dims(A_norm, axis=0) 
# print(A_norm)

# # predict 
# y_norm_pred = model.predict(A_norm)    
# # extractor = keras.Model(inputs=model.inputs, outputs=[layer.output for layer in model.layers])
# # activations = extractor.predict(A_norm)

# # # Save each layer's output to a CSV
# # for i, activation in enumerate(activations):
# #     layer_name = model.layers[i].name
# #     file_name = f"layer_{i}_{layer_name}_output.csv"
# #     file_path = os.path.join(OUTPUT_FOLDER, file_name)

# #     flattened = activation.flatten()
# #     np.savetxt(file_path, flattened, delimiter=",")


# # un-normalize
# invA_pred = (y_norm_pred * (y_std + EPS) + y_mean).reshape(96, 96)

# # print
# np.set_printoptions(precision=15, suppress=True)
# print("Predicted A⁻¹:\n", invA_pred)

# # compare
# true_invA = np.loadtxt(A_PATH.replace('A_','inv_A_'), delimiter=',')
# err = np.linalg.norm(invA_pred - true_invA)
# print(f"Frobenius ‖pred – true‖ = {err:.5e}")



































































#!/usr/bin/env python3
import numpy as np
from tensorflow import keras
import os

MODEL_PATH     = '../MLP_3.keras'
STATS_CSV_PATH = '../MLP_3.csv'
A_PATH         = '../data_200/A_199.csv'
# added paths for b and true x
B_PATH         = '../data_200/b_199.csv'
X_TRUE_PATH    = '../data_200/x_199.csv'
EPS = 1e-8
OUTPUT_FOLDER = "./py_layer_outputs"

# input model
model = keras.models.load_model(MODEL_PATH)

# load csv filess
X_mean = X_std = None
# y_mean = y_std = None
with open(STATS_CSV_PATH, 'r') as f:
    for line in f:
        key, rest = line.split(':', 1)
        vals_str = rest.strip().lstrip('[').rstrip(']')
        arr = np.fromstring(vals_str, sep=',', dtype=np.float32)
        if key.strip() == 'input_mean':
            X_mean = arr
        elif key.strip() == 'input_std':
            X_std = arr
        # elif key.strip() == 'output_mean':
        #     y_mean = arr
        # elif key.strip() == 'output_std':
        #     y_std = arr

# sanity check
for name, v in [('input_mean', X_mean), ('input_std', X_std)]: 
    assert v is not None, f"{name} not found in {STATS_CSV_PATH}"
#    ('output_mean', y_mean), ('output_std', y_std)

# load data
A = np.loadtxt(A_PATH, delimiter=',', dtype=np.float32)
A_flat = A.ravel()

# load b and true x
b       = np.loadtxt(B_PATH,       delimiter=',', dtype=np.float32)
x_true  = np.loadtxt(X_TRUE_PATH,  delimiter=',', dtype=np.float32)

# normalize A
A_norm = (A_flat - X_mean) / (X_std + EPS)
A_norm = np.expand_dims(A_norm, axis=0)
# print(A_norm)

# prepare b input
b_flat  = b.ravel()
b_input = np.expand_dims(b_flat, axis=0)

# predict x
# y_norm_pred = model.predict(A_norm)    
# invA_pred    = (y_norm_pred * (y_std + EPS) + y_mean).reshape(96, 96)
x_pred = model.predict([A_norm, b_input]).flatten()

# extractor = keras.Model(inputs=model.inputs, outputs=[layer.output for layer in model.layers])
# activations = extractor.predict(A_norm)

# # Save each layer's output to a CSV
# for i, activation in enumerate(activations):
#     layer_name = model.layers[i].name
#     file_name = f"layer_{i}_{layer_name}_output.csv"
#     file_path = os.path.join(OUTPUT_FOLDER, file_name)
#
#     flattened = activation.flatten()
#     np.savetxt(file_path, flattened, delimiter=",")

# print outputs
np.set_printoptions(precision=15, suppress=True)
# print("Predicted A⁻¹:\n", invA_pred)
print("predicted x:\n", x_pred)

# compare to true x
# true_invA = np.loadtxt(A_PATH.replace('A_','inv_A_'), delimiter=',')
# err = np.linalg.norm(invA_pred - true_invA)
err = np.linalg.norm(x_pred - x_true)
# print(f"Frobenius ‖pred – true‖ = {err:.5e}")
print(f"frobenius ‖pred – true‖ = {err:.5e}")
