from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Perceptron
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier

from sklearn.metrics import accuracy_score

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd

def plot_decision_regions(x, y, classifier, test_idx=None, resolution=0.02):
  markers = ('s', '1', 'o', '^', 'v','+','p','d')
  colors = ('red','blue','lightgreen','gray','cyan','brown','orange')
  cmap = ListedColormap(colors[:len(np.unique(y))])

  x1_min, x1_max = x[:,0].min() - 1, x[:,0].max() + 1
  x2_min, x2_max = x[:,1].min() - 1, x[:,1].max() + 1
  xx1, xx2 = np.meshgrid(np.arange(x1_min, x1_max, resolution),
                         np.arange(x2_min, x2_max, resolution))
  z = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
  z = z.reshape(xx1.shape)
  plt.contourf(xx1, xx2, z, alpha=0.4, cmap=cmap)
  plt.xlim(xx1.min(), xx1.max())
  plt.ylim(xx2.min(), xx2.max())

  # plot all samples
  for idx, c1 in enumerate(np.unique(y)):
    plt.scatter(x=x[y==c1, 0], y=x[y==c1, 1],
                alpha=0.8, c=cmap(idx), marker=markers[idx], label=c1)
  
  #if test_idx:
  #  x_test, y_test = x[test_idx,:], y[test_idx]
  #  plt.scatter(x=x_test[:, 0], y=x_test[:, 1],
  #              alpha=1.0, linewidths=1, c='white', s=55, marker='x', label='test_set')


#iris = datasets.load_iris()
#x = iris.data[:, [2, 3]]
#y = iris.target
df = pd.read_csv('ionomers.data')
phase_mapping = {'Discrete':0, '':1,'Percolated':2}
df['phase'] = df['phase'].map(phase_mapping)

x = df.iloc[0:60, [0, 1]].values
y = df.iloc[0:60, 2].values
print np.unique(y)
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=0)

sc = StandardScaler()
sc.fit(x_train)

x_train_std = sc.transform(x_train)
x_test_std = sc.transform(x_test)

# Linear models
#ppn = Perceptron(n_iter=40,eta0=0.1, random_state=0)
#ppn.fit(x_train_std, y_train)
#y_pred = ppn.predict(x_test_std)

# for LR, the regularization strength (1/C) is used to prevent overfitting (complexitity of the model)
#   increasing 1/C increases the complexity -- overfitting (tight decision boundaries)
# see page 114 for good explanation
#lr = LogisticRegression(C=1000.0, random_state=0)
#lr.fit(x_train_std, y_train)
#y_pred = lr.predict(x_test_std)

# for SVM, the variable C controls the bias-variation trade-off:
#   increasing C increases the bias and lowers the variance of the model (overfitting, tight decision boundaries)
# linear SVM
#svm = SVC(kernel='linear', C=0.1, random_state=0)

# nonlinear kernel SVM: rbf = Radial Basis Function, gamma control the extent a sample contributes to the model
#svm = SVC(kernel='rbf', C=1.0, gamma=0.8, random_state=0)

#svm.fit(x_train_std, y_train)
#y_pred = svm.predict(x_test_std)

# k-neighbors: nonparametric supervised classification
# p = 2: the Minkowski distance becomes Euclidean distance
#knn = KNeighborsClassifier(n_neighbors=5, p=2, metric='minkowski')
#knn.fit(x_train_std, y_train)
#y_pred = knn.predict(x_test_std)

# decision tree:
tree = DecisionTreeClassifier(criterion='entropy', max_depth=8, random_state=1)
tree.fit(x_train_std, y_train)
y_pred = tree.predict(x_test_std)

print('Misclassified samples: %d' % (y_test != y_pred).sum())
print('Accuracy: %f' % accuracy_score(y_test, y_pred))

x_combined_std = np.vstack((x_train_std, x_test_std))
y_combined = np.hstack((y_train, y_test))
plot_decision_regions(x=x_combined_std, y=y_combined, classifier=tree, test_idx=range(40,55))
plt.xlabel('fq [standardized]')
plt.ylabel('lB [standardized]')
plt.legend(loc='upper left')
plt.show()






