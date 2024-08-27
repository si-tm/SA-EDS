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
from sklearn.metrics import r2_score


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
        # y_data.append(value_dic[(key[0], key[1][1:])]['mean_volume'])
        y_data.append(value_dic[(key[0], key[1][1:])]['sigmoid']["b"]/value_dic[(key[0], key[1][1:])]['sigmoid']["c"])


    x_data = np.array(x_data)
    y_data = np.array(y_data)

    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2)
    return x_train, x_test, y_train, y_test


def bagging_regressor(x_train,x_test,y_train,y_test):
    X = x_train
    y = y_train
    er = BaggingRegressor(estimator=HistGradientBoostingRegressor(),n_estimators=10)

    er = er.fit(X, y)

        # Save the model
    with open(f'/home/user/SA-EDS/saved_model/bagging_model_{type_of_l}_{target}_sigmoid.pkl', 'wb') as f:
        pickle.dump(er, f)
    

    test_predictions = er.predict(x_test)
    train_predictions = er.predict(x_train)
    res = linregress(y_test, test_predictions)

    # if type_of_l == "L1" or type_of_l == "L2":
    #     max_value = 20
    # elif type_of_l == "L3":
    #     max_value = 150

    max_lim = 20
    min_lim = -10

    plt.figure(figsize=(6, 6))
    plt.scatter(y_train, train_predictions, label='Train data')
    plt.scatter(y_test, test_predictions, label='Test data')
    plt.xlabel('True sigmoid value', fontsize=20)
    plt.ylabel('Predicted sigmoid value', fontsize=20)
    plt.axis('equal')
    plt.axis('square')
    plt.xlim([min_lim, max_lim])
    plt.ylim([min_lim, max_lim])
    plt.plot([min_lim, max_lim], [min_lim, max_lim], 'k--', label='Ideal line')

    r2 = r2_score(y_test, test_predictions)
    # Plot the fitted line
    fitted_line = res.intercept + res.slope * np.array([min_lim, max_lim])
    plt.plot([min_lim, max_lim], fitted_line, 'r', label='Fitted line')
    plt.legend()

    # Annotate the R² value on the plot
    plt.text(2, max_lim - 5, f'$R^2 = {r2:.2f}$', fontsize=20, color='red')

    # Display the plot
    plt.show()
    plt.savefig(f'home/user/SA-EDS/results/plot_{type_of_l}_sigmoid.png')  # Save the plot as an image file instead of showing it
    print(f"Plot saved as plot_{type_of_l}_sigmoid.png")


    # Print the fitted line equation and R² score
    print(f"Fitted line: y = {res.intercept} + {res.slope}x")
    print(f"R² score for test data: {r2}")

if __name__ == "__main__": 
    # x_train, x_test, y_train, y_test = get_data_from_f()
    if len(sys.argv) != 3:
        print("usage : python3 predicvVolume.py {L1, L2, L3} target")
    x_train, x_test, y_train, y_test = get_data()
    bagging_regressor(x_train, x_test, y_train, y_test)