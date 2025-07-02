# #####################
# ## 2D CNN VERSION ##
# #####################
# import math
# import numpy as np
# import matplotlib.pyplot as plt
# import scipy as sp
# import tensorflow as tf
# from tensorflow import keras
# import os
# # from mlp import custom_loss

# OUTPUT_FOLDER = "./py_layer_outputs"

# # load model and input matrix
# model = keras.models.load_model('../bin/CNN_1.keras')
# model.summary()
# A = np.loadtxt('../EULER_825/A_824.csv', delimiter=',').astype(np.float32)
# A = A.reshape(-1, 96, 96, 1)

# #####################################
# ## UNCOMMENT FOR PER LAYER OUTPUTS ##
# #####################################
# # extractor = keras.Model(inputs=model.inputs, outputs=[layer.output for layer in model.layers])
# # activations = extractor.predict(A)

# # # Save each layer's output to a CSV
# # for i, activation in enumerate(activations):
# #     layer_name = model.layers[i].name
# #     file_name = f"layer_{i}_{layer_name}_output.csv"
# #     file_path = os.path.join(OUTPUT_FOLDER, file_name)

# #     flattened = activation.flatten()
# #     np.savetxt(file_path, flattened, delimiter=",")

# # perform prediction
# predictions = model.predict(A)
# # print(model.output_shape)
# # np.set_printoptions(threshold=np.inf)
# print(predictions)
































































###################################
## MLP INV_A FLAT ARRAY VERSION ##
###################################
#!/usr/bin/env python3
import numpy as np
from tensorflow import keras
import os
import scipy as sp

MODEL_PATH = "../bin/MLP_LU_1.keras"
STATS_CSV_PATH = "../bin/MLP_LU_1.csv"
A_PATH = "../EULER_825/A_824.csv"
A_INV_PATH = "../EULER_825/A_inv_824.csv"
EPS = 1e-8
OUTPUT_FOLDER = "./py_layer_outputs"

# input model
model = keras.models.load_model(MODEL_PATH)

# load normalization parameters
X_mean = X_std = y_mean = y_std = None
with open(STATS_CSV_PATH, "r") as f:
    for line in f:
        key, rest = line.split(":", 1)
        vals_str = rest.strip().lstrip("[").rstrip("]")
        arr = np.fromstring(vals_str, sep=",", dtype=np.float32)
        if key.strip() == "input_mean":
            X_mean = arr
        elif key.strip() == "input_std":
            X_std = arr
        elif key.strip() == "output_mean":
            y_mean = arr
        elif key.strip() == "output_std":
            y_std = arr
# sanity check
for name, v in [
    ("input_mean", X_mean),
    ("input_std", X_std),
    ("output_mean", y_mean),
    ("output_std", y_std),
]:
    assert v is not None, f"{name} not found in {STATS_CSV_PATH}"

#neural net load and normalize data
X = np.loadtxt(A_PATH, delimiter=",", dtype=np.float32)
inv_A = np.loadtxt(A_INV_PATH, delimiter=",", dtype=np.float32)
X = X.ravel()
X = (X - X_mean) / (X_std + EPS)
X = np.expand_dims(X, axis=0)

# #####################################
# ## UNCOMMENT FOR PER LAYER OUTPUTS ##
# #####################################
# extractor = keras.Model(inputs=model.inputs, outputs=[layer.output for layer in model.layers])
# activations = extractor.predict(X)
# # Save each layer's output to a CSV
# for i, activation in enumerate(activations):
#     layer_name = model.layers[i].name
#     file_name = f"layer_{i}_{layer_name}_output.csv"
#     file_path = os.path.join(OUTPUT_FOLDER, file_name)

#     flattened = activation.flatten()
#     np.savetxt(file_path, flattened, delimiter=",")


# predict and un-normalize
Y = model.predict(X)
Y = (Y * (y_std + EPS) + y_mean)

############
## FOR LU ##
############
#predicted L and U
Y = np.split(Y, 2, axis=1)
Y_L = Y[0].reshape(96, 96)
Y_U = Y[1].reshape(96, 96)
print("\nPredicted L:\n", Y_L)
print("\nPredicted U:\n", Y_U)
#actual L and U
L, U = sp.linalg.lu(inv_A, permute_l=True)
# P, L, U = sp.linalg.lu(A, permute_l=False)
print("\nActual L:\n", L)
print("\nActual U:\n", U)

# #print
# np.set_printoptions(precision=15, suppress=True)
# print("Predicted A⁻¹:\n", Y)































































#####################################
## MLP X WITH B FLAT ARRAY VERSION ##
#####################################
# #!/usr/bin/env python3
# import numpy as np
# from tensorflow import keras
# import os

# MODEL_PATH     = '../MLP_x_1.keras'
# STATS_CSV_PATH = '../MLP_x_1.csv'
# A_PATH         = '../EULER_825/A_0.csv'
# # added paths for b and true x
# # B_PATH         = '../data_200/b_199.csv'
# # X_TRUE_PATH    = '../data_200/x_199.csv'
# EPS = 1e-8
# OUTPUT_FOLDER = "./py_layer_outputs"

# # input model
# model = keras.models.load_model(MODEL_PATH, compile=False)
# model.summary()

# # load csv filess
# X_mean = X_std = None
# y_mean = y_std = None
# with open(STATS_CSV_PATH, 'r') as f:
#     for line in f:
#         key, rest = line.split(':', 1)
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
# for name, v in [('input_mean', X_mean), ('input_std', X_std), ('output_mean', y_mean), ('output_std', y_std)]:
#     assert v is not None, f"{name} not found in {STATS_CSV_PATH}"

# # load data
# A = np.loadtxt(A_PATH, delimiter=',', dtype=np.float32)
# A_flat = A.ravel()

# # load b and true x
# # b       = np.loadtxt(B_PATH,       delimiter=',', dtype=np.float32)
# # x_true  = np.loadtxt(X_TRUE_PATH,  delimiter=',', dtype=np.float32)

# # normalize A
# A_norm = (A_flat - X_mean) / (X_std)
# A_norm = np.expand_dims(A_norm, axis=0)
# # print(A_norm)

# # prepare b input
# # b_flat  = b.ravel()
# # b_input = np.expand_dims(b_flat, axis=0)

# # predict x
# y_norm_pred = model.predict(A_norm)
# invA_pred    = (y_norm_pred * (y_std) + y_mean).reshape(96, 96)
# # x_pred = model.predict([A_norm, b_input]).flatten()

# # # per layer output
# # extractor = keras.Model(inputs=model.inputs, outputs=[layer.output for layer in model.layers])
# # activations = extractor.predict(A_norm)
# # for i, activation in enumerate(activations):
# #     layer_name = model.layers[i].name
# #     file_name = f"layer_{i}_{layer_name}_output.csv"
# #     file_path = os.path.join(OUTPUT_FOLDER, file_name)

# #     flattened = activation.flatten()
# #     np.savetxt(file_path, flattened, delimiter=",")

# # print outputs
# np.set_printoptions(precision=15, suppress=True)
# print("Predicted A⁻¹:\n", invA_pred)
# # print("predicted x:\n", x_pred)

# # compare to true x
# # true_invA = np.loadtxt(A_PATH.replace('A_','inv_A_'), delimiter=',')
# # err = np.linalg.norm(invA_pred - true_invA)
# # err = np.linalg.norm(x_pred - x_true)
# # print(f"Frobenius ‖pred – true‖ = {err:.5e}")
