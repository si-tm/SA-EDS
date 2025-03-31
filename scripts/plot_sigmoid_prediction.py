import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import sys
from classify_seq import make_input_seq as mis
from sklearn.model_selection import train_test_split

def load_and_plot_saved_model(type_of_l, target, x_test, y_test, x_train, y_train):
    # 保存されたモデルをロード
    model_path = f'/home/user/SA-EDS/saved_model/bagging_model_{type_of_l}_{target}_sigmoid.pkl'
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    
    # テストデータで予測
    test_predictions = model.predict(x_test)
    train_predictions = model.predict(x_train)
    res = linregress(y_test, test_predictions)
    
    # 軸の最大値設定
    if type_of_l in ["L1", "L2"]:
        max_value = 40
    elif type_of_l == "L3":
        max_value = 150
    else:
        max_value = max(y_test) * 1.2
    
    # グラフの描画
    plt.figure(figsize=(6,6))
    
    # トレーニングデータを青色でプロット
    plt.scatter(y_train, train_predictions, label="Train Data", alpha=0.6, color='b')
    
    # テストデータを緑色でプロット
    plt.scatter(y_test, test_predictions, label="Test Data", alpha=0.6, color='r')
    
    # 理想的な直線
    plt.plot([0, max_value], [0, max_value], 'g--', label="Ideal Line")
    plt.plot([0, max_value], res.intercept + res.slope * np.array([0, max_value]), 'r', label="Fitted Line")
    
    r_squared = model.score(x_test, y_test)
    plt.title(f"{type_of_l} Volume Prediction (Reloaded Model)\nR²: {r_squared:.2f}")
    plt.xlabel("Actual Volume")
    plt.ylabel("Predicted Volume")
    plt.legend()
    
    # 画像を保存
    save_path = f'/home/user/SA-EDS/results/reloaded_plot_{type_of_l}_{target}_sub_sigmoid.png'
    plt.savefig(save_path)
    print(f"Plot saved as {save_path}, R²: {r_squared:.2f}")
    plt.show()

def get_data(target_lst, type_of_l):
    domain_lst = mis.seq_lst(f"home/user/SA-EDS/conf/input_seq_{type_of_l}.csv")
    domain_seq_dic = {}
    value_dic = {}
    target_lst = [x for x in target_lst if x]

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


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python3 reloadModelPlot.py {L1, L2, L3} target1 target2 ...")
        sys.exit(1)
    
    type_of_l = sys.argv[1]
    target_lst = sys.argv[2:]
    target_lst = [x for x in target_lst if x]
    target_name = "-".join(target_lst)
    
    x_train, x_test, y_train, y_test = get_data(target_lst, type_of_l)
    
    load_and_plot_saved_model(type_of_l, target_name, x_test, y_test, x_train, y_train)
