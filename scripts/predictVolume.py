import tensorflow as tf
print("TensorFlow version:", tf.__version__)

from tensorflow.keras.layers import Dense
from tensorflow.keras import Model
from tensorflow.keras import layers
import sys
import matplotlib.pyplot as plt
import pickle
from scipy.stats import linregress
from sklearn.model_selection import train_test_split
from datetime import datetime
import numpy as np
from sklearn.ensemble import BaggingRegressor
from sklearn.datasets import make_regression

sys.path.append("../measuring_volume")
sys.path.append("../common")
import pickle
from classify_seq import make_input_seq as mis


from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor, GradientBoostingRegressor, StackingRegressor, VotingRegressor
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
import numpy as np
import tensorflow as tf

# type_of_l = "L1"
type_of_l = sys.argv[1]
target = sys.argv[2]

def get_data():
    # [0, 1, ... 0] + temperature
    f1 = open(f"home/user/SA-EDS/dataset/x_{target}_{type_of_l}.pkl", "rb")
    f2 = open(f"home/user/SA-EDS/dataset/{type_of_l}_data_{target}.pkl", "rb")
    domain_seq_dic = pickle.load(f1)
    value_dic = pickle.load(f2)
    domain_lst = mis.seq_lst(f"home/user/SA-EDS/conf/input_seq_{type_of_l}.csv")

    f1.close()
    f2.close()

    x_data = []
    y_data = []

    for key in domain_seq_dic:
        new_x = []
        for domain in domain_lst:
            new_x.append(domain_seq_dic[key]["domain"][domain])
        new_x.append(int(key[0]))
        new_x.append(float(domain_seq_dic[key]['eigenValue_2']))
        x_data.append(new_x)
        y_data.append(value_dic[(key[0], key[1][1:])]['mean_volume'])

    x_data = np.array(x_data)
    y_data = np.array(y_data)

    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2)
    return x_train, x_test, y_train, y_test

# def get_data_from_f(f_x2=["dataset/x_QD_1_L1.pkl"], f_y2=["dataset/L1_data_QD_1.pkl"]):
#     # [0, 1, ... 0] + temperature
#     f_x1 = open(f"dataset/x_initial_{type_of_l}.pkl", "rb")
#     f_y1= open(f"dataset/{type_of_l}_data_initial.pkl", "rb")
#     domain_seq_dic = pickle.load(f_x1)
#     value_dic = pickle.load(f_y1)
#     for fx in f_x2:
#         fx = open(fx, "rb")
#         domain_seq_dic.update(pickle.load(fx))
#     for fy in f_y2:
#         fy = open(fy, "rb")
#         value_dic.update(pickle.load(fy))
#     domain_lst = mis.seq_lst(f"conf/input_seq_{type_of_l}.csv")

#     f_x1.close()
#     f_y1.close()

#     x_data = []
#     y_data = []

#     for key in domain_seq_dic:
#         new_x = []
#         print(domain_seq_dic[key])
#         for domain in domain_lst:
#             new_x.append(domain_seq_dic[key]["domain"][domain])
#         new_x.append(int(key[0]))
#         new_x.append(float(domain_seq_dic[key]['eigenValue_2']))
#         x_data.append(new_x)
#         y_data.append(value_dic[key]['mean_volume'])

#     x_data = np.array(x_data)
#     y_data = np.array(y_data)

#     x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2)
#     return x_train, x_test, y_train, y_test



def bagging_regressor(x_train,x_test,y_train,y_test):
    X = x_train
    y = y_train
    er = BaggingRegressor(estimator=RandomForestRegressor(n_estimators=10),n_estimators=10)
    er = er.fit(X, y)

        # Save the model
    with open(f'/home/user/SA-EDS/saved_model/bagging_model_{type_of_l}_{target}.pkl', 'wb') as f:
        pickle.dump(er, f)
    

    test_predictions = er.predict(x_test)
    train_predictions = er.predict(x_train)
    res = linregress(y_test, test_predictions)

    if type_of_l == "L1" or type_of_l == "L2":
        max_value = 20
    elif type_of_l == "L3":
        max_value = 150

    plt.figure(figsize=(6,6))
    plt.scatter(y_train, train_predictions)
    plt.scatter(y_test, test_predictions)
    plt.xlabel(f'{type_of_l} True average of DNA volume')
    plt.ylabel(f'{type_of_l} Prediction of DNA volume average')
    plt.axis('equal')
    plt.axis('square')
    plt.xlim([0,max_value])
    plt.ylim([0,max_value])
    _ = plt.plot([0,max_value], [0,max_value])

    print(er.score(x_test, y_test))

    plt.plot([0,max_value], res.intercept + res.slope*np.array([0, max_value]), 'r', label='fitted line')

    plt.show()
    plt.savefig(f'home/user/SA-EDS/results/plot_{type_of_l}_{target}.png')  # Save the plot as an image file instead of showing it
    print(f"Plot saved as plot_{type_of_l}_{target}.png")

if __name__ == "__main__": 
    # x_train, x_test, y_train, y_test = get_data_from_f()
    if len(sys.argv) != 3:
        print("usage : python3 predicvVolume.py {L1, L2, L3} target")
    x_train, x_test, y_train, y_test = get_data()
    bagging_regressor(x_train, x_test, y_train, y_test)