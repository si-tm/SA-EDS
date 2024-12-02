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
import numpy as np
from sklearn.ensemble import BaggingRegressor, RandomForestRegressor
import sys

sys.path.append("../measuring_volume")
sys.path.append("../common")
from classify_seq import make_input_seq as mis

type_of_l = sys.argv[1]

def get_data(target_lst):
    domain_lst = mis.seq_lst(f"home/user/SA-EDS/conf/input_seq_{type_of_l}.csv")
    domain_seq_dic = {}
    value_dic = {}

    for target in target_lst:
        f1 = open(f"home/user/SA-EDS/dataset/x_{target}_{type_of_l}.pkl", "rb")
        f2 = open(f"home/user/SA-EDS/dataset/{type_of_l}_data_{target}.pkl", "rb")
        domain_seq_dic.update(pickle.load(f1))
        value_dic.update(pickle.load(f2))

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

    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.1)
    return x_train, x_test, y_train, y_test

def bagging_regressor(x_train, x_test, y_train, y_test, target):
    er = BaggingRegressor(estimator=RandomForestRegressor(n_estimators=10), n_estimators=10)
    er.fit(x_train, y_train)

    with open(f'/home/user/SA-EDS/saved_model/bagging_model_{type_of_l}_{target}.pkl', 'wb') as f:
        pickle.dump(er, f)
    test_predictions = er.predict(x_test)
    train_predictions = er.predict(x_train)
    res = linregress(y_test, test_predictions)

    if type_of_l == "L1" or type_of_l == "L2":
        max_value = 20
        max_value = 40
    elif type_of_l == "L3":
        max_value = 150

    plt.figure(figsize=(6,6))
    plt.scatter(y_train, train_predictions, label="Train Data", alpha=0.6)
    plt.scatter(y_test, test_predictions, label="Test Data", alpha=0.6)
    plt.plot([0, max_value], [0, max_value], 'g--', label="Ideal Line")  # Ideal line in green dashed
    plt.plot([0, max_value], res.intercept + res.slope * np.array([0, max_value]), 'r', label="Fitted Line")  # Fitted line in red
    plt.xlabel(f'{type_of_l} True average of DNA volume')
    plt.ylabel(f'{type_of_l} Prediction of DNA volume average')
    plt.xlim([0, max_value])
    plt.ylim([0, max_value])
    plt.legend(loc="upper left")

    r_squared = er.score(x_test, y_test)
    plt.title(f"{type_of_l} Volume Prediction\nR²: {r_squared:.2f}")

    plt.savefig(f'home/user/SA-EDS/results/plot_{type_of_l}_{target}.png')  # Save plot as image
    print(f"Plot saved as plot_{type_of_l}_{target}.png, R²: {r_squared:.2f}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage : python3 predictVolume.py {L1, L2, L3} target1 target2 ...")
    
    target_lst = []

    print(sys.argv)
    
    for i, argv in enumerate(sys.argv):
        if i > 1:
            target_lst.append(argv)
    target_lst = list(filter(bool, target_lst))
    print(target_lst)

    x_train, x_test, y_train, y_test = get_data(target_lst)
    bagging_regressor(x_train, x_test, y_train, y_test, ("-").join(target_lst))
