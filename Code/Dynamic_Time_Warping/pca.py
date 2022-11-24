import os
import sys
import csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn import svm
from sklearn.preprocessing import normalize as norm
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import scale
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

rootdir = "../PCA_Data/"

def main(mean, plot):

    #
    #   Currently, reduces down the number of phonemes groups - make into func and add another to reduce down
    #   coefficients should be ["Coefficient
    #
    feats = ["Coefficient " + str(i) for i in range(0, 14)] + ["Phoneme"]
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


    for subdir, dirs, files in os.walk(rootdir):
        for tmp in files:
            if tmp.endswith('.txt'):
                with open(subdir + "/" + tmp, 'r') as file:
                    final_data = pd.read_csv(file, names=feats, delimiter=",")

    # features = ["Coefficient " + str(i) for i in range(0, 24)]
    feats = ["Coefficient " + str(i) for i in range(0, 14)]
    fig = px.scatter_matrix(
        final_data,
        dimensions=feats,
        color="Phoneme"
    )
    fig.update_traces(diagonal_visible=False)
    if plot == ['0']:
        fig.show() # First plot of all phonemes
  
    print("Final data :: \n")
    print(final_data[feats])
    # PCA
    print("PCA")
    pca = PCA(n_components=2)
    components = pca.fit_transform(final_data[feats])

    # K-means
    print("Kmeans");
    kmeans = KMeans(n_clusters=38, random_state=0).fit(components)
    k_labels = list(set(final_data['Phoneme']))

    # SVM
    print("SVM")
    print("Normalising")
    # normed = norm(final_data[feats])
    scaler = MinMaxScaler()
    normed = scaler.fit_transform(final_data[feats])
    print("Fitting")
    # clf = OneVsRestClassifier(SVC(kernel='rbf', verbose=True, probability=True, class_weight='auto'), n_jobs=-1).fit(normed, final_data['Phoneme'])
    clf = svm.LinearSVC(verbose=1, max_iter=10000)
    clf.fit(normed, final_data['Phoneme'])

    print("Weights")
    print(clf.coef_)
    with open(rootdir + "/res/weights.org", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for i in clf.coef_:
            writer.writerow(i)
    print("Constants")
    print(clf.intercept_)
    with open(rootdir + "/res/constants.org", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')       
        writer.writerow(clf.intercept_)
    

    
    if mean == ['0']:
        fig = px.scatter(kmeans.cluster_centers_, x=0, y=1, color=k_labels)
    else:
        fig = px.scatter(components, x=0, y=1, color=final_data["Phoneme"])
        
    if plot == ['0']:
        fig.show() # Second plot
    
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    fig = px.scatter(components, x=0, y=0, color=final_data["Phoneme"])
    print(final_data['Phoneme'])
    print(loadings)
    for i, feature in enumerate(feats):
        fig.add_shape(
            type='line',
            x0=0, y0=0,
            x1=loadings[i, 0],
            y1=loadings[i, 1]
        )
        fig.add_annotation(
            x=loadings[i, 0],
            y=loadings[i, 1],
            ax=0, ay=0,
            xanchor="center",
            yanchor="bottom",
            text=feature,
        )
    if plot == ['0']:
        fig.show() # Plot with the lines

    PC_values = np.arange(pca.n_components_) + 1
    plt.plot(PC_values, pca.explained_variance_ratio_, 'o-', linewidth=2, color='blue')
    plt.title('Scree Plot')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    if plot == ['0']:
        plt.show() # Scree plot

    if not os.path.exists(rootdir + "/res/"):
        os.makedirs(rootdir + "/res/")
    if not os.path.exists(rootdir + "/fig/"):
        os.makedirs(rootdir + "/fig/")

    # for i in range(0, len(final_data)):
    #     set = final_data[i][0:6]
    #     for j in range(0, len(final_data[i])):
    #         final_data[i][j] = loadings[j] * set[i][j]
    #
    # data_set = final_data
    # print(data_set)

    # with open(rootdir + "/res/final.pca", 'w', newline='') as csvfile:
    #     writer = csv.writer(csvfile, delimiter=',')
    #     for i in data_set:
    #         writer.writerow(i)
    # plt.savefig("/figs/_pca.png")
    # plt.clf()
    with open(rootdir + "/res/loadings.org", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        for i in loadings:
            writer.writerow(i)
    with open(rootdir + "/res/kmeans.org", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        for i in kmeans.cluster_centers_:
            writer.writerow(i)
        


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
