import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def open_file():
    with open('confusion_matrix.txt', 'r') as file:
        data = file.read().split("\n")
        labels = list()
        for i in range(0, len(data)):
            tmp = data[i].split("|")
            if tmp[0] == '':
                continue
            labels.append(tmp[0])
            tmp.pop(0)
            data[i] = tmp
            data[i].pop(len(data[i]) - 1)
    return data, labels


def gen_matrix():
    data, labels = open_file()
    data.pop(len(data) - 1)
    data = np.array(data, dtype=float)
    ax = plt.subplot()
    sns.heatmap(data, xticklabels=labels, yticklabels=labels)
    # sns.heatmap(data, annot=True, ax=ax)
    # sns.heatmap(data, linewidths=0.5, linecolor="red")
    plt.show()


if __name__ == "__main__":
    gen_matrix()
