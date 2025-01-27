
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import pickle
import numpy as np

"""
        QD version                                     feature-2
0  preliminary L1          $\frac{(\frac{b}{c} + 3)}{(5 + 3)}$
1   seconderly L1         $\frac{(\frac{b}{c} - 5)}{(23 - 5)}$
2  preliminary L2          $\frac{(\frac{b}{c} + 3)}{(5 + 3)}$
3   seconderly L2         $\frac{(\frac{b}{c} - 5)}{(15 - 5)}$
4  preliminary L3          $\frac{(\frac{b}{c} + 3)}{(5 + 3)}$
5   seconderly L3        $\frac{(\frac{b}{c} - 40)}{(100 - 40)}$
"""

"""
x = (b/c + 3)/(5 + 3)
x*(5 + 3) = b/c + 3
x*(5 + 3) - 3 = b/c
"""

def denormalize_bc(bc, type_of_l, target):
    denormalize_val = {
        "initial": {
            "L1": [5, 3],
            "L2": [5, 3],
            "L3": [5, 3]
        },
        "second": {
            "L1": [23, 5],
            "L2": [15, 5],
            "L3": [100, 40],
        },
    }
    print(denormalize_val[target][type_of_l][1], denormalize_val[target][type_of_l][0], denormalize_val[target][type_of_l][1])
    return bc*(denormalize_val[target][type_of_l][0] - denormalize_val[target][type_of_l][1]) + denormalize_val[target][type_of_l][1]

def comp_stability_in_QD(type_of_l="L1"):
    base_dir = "/home/user/SA-EDS/"
    # Load the data
    f_initial = open(f'{base_dir}results/qd_output_int_initial_{type_of_l}.csv')
    f_second = open(f'{base_dir}results/qd_output_int_second_{type_of_l}.csv')

    initial_bc = []
    second_bc = []

    for i, l in enumerate(f_initial):
        if i == 0:
            continue
        bc = l.split(",")[3]
        bc = denormalize_bc(float(bc), type_of_l, "initial")
        initial_bc.append(bc)
    for i, l in enumerate(f_second):
        if i == 0:
            continue
        bc = l.split(",")[3]
        bc = denormalize_bc(float(bc), type_of_l, "second")
        second_bc.append(bc)

    # Convert the strings to floats
    initial_bc = [float(x) for x in initial_bc]
    second_bc = [float(x) for x in second_bc]

    # Convert the lists to Series
    initial_bc = pd.Series(initial_bc)
    second_bc = pd.Series(second_bc)

    # Create a DataFrame for plotting
    data = pd.DataFrame({
        'Initial': initial_bc,
        'Second': second_bc
    })

    # Plot the boxplot
    plt.figure(figsize=(10, 6))
    data.boxplot()
    plt.title('Boxplot of b/c values')
    plt.ylabel('b/c')
    plt.show()
    plt.savefig(f'{base_dir}results/comp_stability_{type_of_l}.png')

    return data

def main():
    if len(sys.argv) != 4:
        print("usage: python3 comp_features.py <type of l>")
    type_of_l = sys.argv[1]
    comp_stability(type_of_l)

def plot_all_in_QD():
    data1 = comp_stability_in_QD("L1")
    data2 = comp_stability_in_QD("L2")
    data3 = comp_stability_in_QD("L3")

    # Concatenate the dataframes correctly
    data = pd.DataFrame({
        'Initial L1': data1['Initial'],
        'Second L1': data1['Second'],
        'Initial L2': data2['Initial'],
        'Second L2': data2['Second'],
        'Initial L3': data3['Initial'],
        'Second L3': data3['Second']
    })

    # plot the boxplot of all the data in each type of L
    # all fontsize is 20
    
    # plt.legend(fontsize=35)
    plt.figure(figsize=(12, 8))
    data.boxplot()
    plt.title('stability b/c plotted in QD', fontsize=35)
    plt.ylabel('b/c', fontsize=35)
    plt.show()
    base_dir = "/home/user/SA-EDS/"
    plt.savefig(f'{base_dir}results/comp_stability_in_QD.png')


def comp_stability(type_of_l="L1"):
    base_dir = "/home/user/SA-EDS/"
    # Load the data
    f_initial = open(f"{base_dir}dataset/{type_of_l}_data_int_initial.pkl", "rb")
    f_second = open(f"{base_dir}dataset/{type_of_l}_data_int_second_each.pkl", "rb")

    initial_bc = []
    second_bc = []

    initial_data = pickle.load(f_initial)
    second_data = pickle.load(f_second)

    for i, l in enumerate(initial_data):
        if "sigmoid" in initial_data[l]:
            if np.any(initial_data[l]["sigmoid"]["c"] != 0):
                initial_bc.append(initial_data[l]["sigmoid"]["b"]/initial_data[l]["sigmoid"]["c"])

    for i, l in enumerate(second_data):
        if "sigmoid" in second_data[l]:
            if np.any(second_data[l]["sigmoid"]["c"] != 0):
                print(second_data[l]["sigmoid"])
                try:
                    second_bc.append(second_data[l]["sigmoid"]["b"]/second_data[l]["sigmoid"]["c"][2])
                except:
                    second_bc.append(second_data[l]["sigmoid"]["b"]/second_data[l]["sigmoid"]["c"])

                # print(second_data[l]["sigmoid"]["b"], second_data[l]["sigmoid"]["b"]/second_data[l]["sigmoid"]["c"])

    # Convert the lists to Series
    initial_bc = pd.Series(initial_bc)
    second_bc = pd.Series(second_bc)
    print(second_bc)
    # Create a DataFrame for plotting
    data = pd.DataFrame({
        'Initial': initial_bc,
        'Second': second_bc
    })

    # print(data)
    
    return data


def plot_all():
    data1 = comp_stability("L1")
    data2 = comp_stability("L2")
    data3 = comp_stability("L3")

    # print(data1)
    # print(data2)
    # print(data3)

    # Create a DataFrame for plotting
    data = pd.DataFrame({
        'Initial L1': data1['Initial'],
        'Second L1': data1['Second'],
        'Initial L2': data2['Initial'],
        'Second L2': data2['Second'],
        'Initial L3': data3['Initial'],
        'Second L3': data3['Second']
    })

    # 

    print("data second l3: ", data['Second L3'])
    # print("data initial l2: ", data['Initial L2'])
    # data_cleaned = data.dropna()

    # plot the boxplot of all the data in each type of L
    # all fontsize is 20
    plt.legend(fontsize=20)
    data.boxplot()
    plt.title('stability b/c of QD result', fontsize=20)
    plt.ylabel('b/c', fontsize=20)
    plt.show()
    base_dir = "/home/user/SA-EDS/"
    plt.savefig(f'{base_dir}results/comp_stability_oxDNA.png')

if __name__ == '__main__':
    # main()
    plot_all_in_QD()
    # plot_all()

