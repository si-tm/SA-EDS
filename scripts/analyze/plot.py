import pickle
import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
    type_of_l = "L1"
    target = "int_initial"

    type_of_l = sys.argv[1]
    target = sys.argv[2]
    base_dir = "/home/user/SA-EDS/"
    volume_dir = base_dir + f"dataset/{type_of_l}_data_{target}.pkl"
    eigen_dir = base_dir + f"dataset/x_{target}_{type_of_l}.pkl"
    
    # Load data
    with open(volume_dir, "rb") as f:
        volume_data = pickle.load(f)
    with open(eigen_dir, "rb") as f:
        eigen_data = pickle.load(f)

    volume_lst = []
    eigen_lst = []

    # Extract data
    for i, l in enumerate(volume_data):
        if l[0] == "277":
            volume_lst.append(volume_data[l]["mean_volume"])
            l = list(l)  # Convert tuple to list
            l[1] = "/" + l[1]  # Modify the second element
            l = tuple(l)  # Optionally, convert back to a tuple if needed
            eigen_lst.append(eigen_data[l]['eigenValue_2'])

    # Calculate Pearson correlation coefficient
    correlation = np.corrcoef(volume_lst, eigen_lst)[0, 1]
    print(f"Pearson correlation coefficient: {correlation}")

    # Plot the data
    plt.figure(figsize=(8, 6))
    plt.plot(volume_lst, eigen_lst, marker='o', linestyle='', color='b')  # Plot points only

    # Fit a linear regression line
    coeffs = np.polyfit(volume_lst, eigen_lst, 1)
    poly_eqn = np.poly1d(coeffs)
    plt.plot(volume_lst, poly_eqn(volume_lst), color='r', linestyle='--', label="Regression Line")

    # Add labels and grid
    plt.title("Volume vs. Eigenvalue 2")
    plt.xlabel("Mean Volume")
    plt.ylabel("Eigenvalue 2")
    plt.grid(True)
    plt.legend()

    # Save the plot
    plt.savefig(base_dir + f"fig/volume_vs_eigen_with_trend_{type_of_l}_{target}.png")
    plt.show()

if __name__ == '__main__':
    main()
