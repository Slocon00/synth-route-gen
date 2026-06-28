from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import warnings


warnings.filterwarnings('ignore', category=RuntimeWarning)

class XMeans:
    '''
    I removed useless parameters (such as max_iter, tol,...) and also I removed distance_function as a parameter to
    XMeans and consequently to KMeans, because it is removed from the latest versions of SKlearn.

    '''
    def __init__(self, kmin, kmax):
        self.kmin = kmin
        self.kmax = kmax
        self.init = 'k-means++'

        self.k_ = -1
        self.cluster_centers_ = None
        self.labels_ = None
        self.inertia_ = None


    def fit(self, data):
        #print("I am using new XMeans written by Omid")
        cluster_centers = []
        k = self.kmin
        #print("** ", self.kmin, self.kmax)
        while k <= self.kmax:
            kmeans = self._fit(k, data, cluster_centers)

            centroids = kmeans.cluster_centers_

            centroid_distances = euclidean_distances(centroids)
            centroid_distances += np.diag([np.inf] * k)
            min_centroid_distances = centroid_distances.min(axis = -1)

            labels = kmeans.labels_

            cluster_centers = []
            for i, centroid in enumerate(centroids):
                direction = np.random.random(centroid.shape)

                vector = direction / np.sqrt(np.dot(direction,direction)) * min_centroid_distances[i]

                new_point1 = centroid + vector
                new_point2 = centroid - vector

                label_index = (labels == i)
                points = data[label_index]

                if len(np.unique(points, axis=0)) == 1:
                    cluster_centers.append(centroid)
                    continue

                new_kmeans = self._fit(2, points, np.asarray([new_point1, new_point2]))
                new_labels = new_kmeans.labels_
                cluster1 = points[new_labels == 0]
                cluster2 = points[new_labels == 1]

                if len(cluster1) == 0 or len(cluster2) == 0:
                    cluster_centers.append(centroid)
                    continue

                bic_parent = XMeans.bic([points], [centroid])
                bic_child = XMeans.bic([cluster1, cluster2], new_kmeans.cluster_centers_)

                if bic_child > bic_parent:
                    cluster_centers.extend(new_kmeans.cluster_centers_)
                else:
                    cluster_centers.append(centroid)
            if k==len(cluster_centers):
                break
            k = len(cluster_centers)

        #print('count of iteration: {}'.format(k))

        # Final K-Means with the best k.
        if len(cluster_centers) == 0:
            cluster_centers = self.init
        else:
            cluster_centers = np.asarray(cluster_centers)
        result = KMeans(k, init=cluster_centers, verbose=False).fit(data)

        self.cluster_centers_ = result.cluster_centers_
        self.labels_ = result.labels_
        self.inertia_ = result.inertia_
        self.k_ = k
        #print("** ", self.k_, self.inertia_)

        return self




    def _fit(self, k, data, centroids):
        if len(centroids) == 0:
            centroids = self.init
        else:
            centroids = np.asarray(centroids)
        result = KMeans(k, init=centroids, verbose=False).fit(data)
        #print(result)
        return result

    @classmethod
    def bic(cls, clusters, centroids):
        R = sum([len(cluster) for cluster in clusters])
        M = clusters[0][0].shape[0]
        K = len(centroids)
        log_likelihood = XMeans._log_likelihood(R, M, clusters, centroids)
        num_params = XMeans._free_params(K, M)
        return log_likelihood - num_params / 2.0 * np.log(R)

    @classmethod
    def _log_likelihood(cls, R, M, clusters, centroids):
        ll = 0
        var = XMeans._variance(R, M, clusters, centroids)
        #print('estimate {}'.format(var))
        for cluster in clusters:
            R_n = len(cluster)
            t1 = R_n * np.log(R_n)
            t2 = R_n * np.log(R)
            t3 = R_n * M / 2.0 * np.log(2.0 * np.pi * var)
            t4 = M * (R_n-1.0) / 2.0
            ll += t1 - t2 - t3 - t4
        return ll

    @classmethod
    def _variance(cls, R, M, clusters, centroids):
        K = len(centroids)
        denom = float((R - K) * M)
        s = 0
        for cluster, centroid in zip(clusters, centroids):
            distances = euclidean_distances(cluster, [centroid])
            s += (distances * distances).sum()
        return s / denom

    @classmethod
    def _free_params(cls, K, M):
        return K * (M+1)