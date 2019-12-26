from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, cdist, squareform
from random import sample 
import pdb
np.random.seed(6)

class GPR():
   
    def __init__(self, kernel, optimizer='L-BFGS-B', noise_var=1e-8):
        self.kernel = kernel
        self.noise_var = noise_var
        self.optimizer = optimizer
    
    
    # 'Public' methods
    def sample_prior(self, X_test, n_samples=1):
        y_mean = np.zeros(X_test.shape[0])
        y_cov = self.kernel(X_test)
        return self._sample_multivariate_gaussian(y_mean, y_cov, n_samples)
    
    
    def sample_posterior(self, X_train, y_train, X_test, n_samples=1):
        # compute alpha
        K = self.kernel(X_train)
        K[np.diag_indices_from(K)] += self.noise_var
        L = self._cholesky_factorise(K)
        alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train))

        # Compute posterior mean
        K_trans = self.kernel(X_test, X_train)
        y_mean = K_trans.dot(alpha)

        # Compute posterior covariance
        v = np.linalg.solve(L, K_trans.T)  # L.T * K_inv * K_trans.T
        y_cov = self.kernel(X_test) - np.dot(v.T, v)

        return self._sample_multivariate_gaussian(y_mean, y_cov, n_samples), y_mean, y_cov
    
    
    def log_marginal_likelihood(self, X_train, y_train, theta, noise_var=None):
    
        if noise_var is None:
            noise_var = self.noise_var

        # Build K(X, X)
        self.kernel.theta = theta
        K = self.kernel(X_train)    
        K[np.diag_indices_from(K)] += noise_var

        # Compute L and alpha for this K (theta).
        L = self._cholesky_factorise(K)
        alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train))

        # Compute log marginal likelihood.
        log_likelihood = -0.5 * np.dot(y_train.T, alpha)
        log_likelihood -= np.log(np.diag(L)).sum()
        log_likelihood -= K.shape[0] / 2 * np.log(2 * np.pi)

        return log_likelihood
    
    
    def optimize(self, X_train, y_train):
    
        def obj_func(theta, X_train, y_train):
                return -self.log_marginal_likelihood(X_train, y_train, theta)

        results = minimize(obj_func, 
                           self.kernel.theta, 
                           args=(X_train, y_train), 
                           method=self.optimizer, 
                           jac=None,
                           bounds=self.kernel.bounds)

        # Store results of optimization.
        self.max_log_marginal_likelihood_value = -results['fun']
        self.kernel.theta_MAP = results['x']

        return results['success']
    
    
    # 'Private' helper methods
    def _cholesky_factorise(self, y_cov):
        try:
            L = np.linalg.cholesky(y_cov)
        except np.linalg.LinAlgError as e:
            e.args = ("The kernel, %s, is not returning a" 
                      "positive definite matrix. Try increasing"
                      " the noise variance of the GP or using"
                      " a larger value for epsilon. "
                      % self.kernel,) + e.args
            raise
        return L
    
    
    def _sample_multivariate_gaussian(self, y_mean, y_cov, n_samples=1, epsilon=1e-10):
        y_cov[np.diag_indices_from(y_cov)] += epsilon  # for numerical stability
        L = self._cholesky_factorise(y_cov)
        u = np.random.randn(y_mean.shape[0], n_samples)
        z = np.dot(L, u) + y_mean[:, np.newaxis]
        return z



class Linear():
    def __init__(self, signal_variance=1.0, signal_variance_bounds=(1e-5, 1e5)):
        self.theta = [signal_variance]
        self.bounds = [signal_variance_bounds]
    def __call__(self, X1, X2=None):
        if X2 is None:
            K = self.theta[0] * np.dot(X1, X1.T)
        else:
            K = self.theta[0] * np.dot(X1, X2.T)
        return K
    

class SquaredExponential():
    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5)):
        self.theta = [length_scale]
        self.bounds = [length_scale_bounds]
    def __call__(self, X1, X2=None):
        if X2 is None:
            # K(X1, X1) is symmetric so avoid redundant computation using pdist.
            dists = pdist(X1 / self.theta[0], metric='sqeuclidean')
            K = np.exp(-0.5 * dists)
            K = squareform(K)
            np.fill_diagonal(K, 1)
        else:
            dists = cdist(X1 / self.theta[0], X2 / self.theta[0], metric='sqeuclidean')
            K = np.exp(-0.5 * dists)
        return K


class RationalQuadratic():
    def __init__(self, length_scale=1.0, alpha=1, length_scale_bounds=(1e-5, 1e5), alpha_bounds=(1e-5,1e5)):
        self.theta = [length_scale, alpha]
        self.bounds = [length_scale_bounds, alpha_bounds]
    def __call__(self, X1, X2=None):
        if X2 is None:
            # K(X1, X1) is symmetric so avoid redundant computation using pdist.
            dists = pdist(X1 / self.theta[0]*np.sqrt(self.theta[1]), metric='sqeuclidean')
            K = (1 + dists)**(-self.theta[1])
            K = squareform(K)
            np.fill_diagonal(K, 1)
        else:
            dists = cdist(X1 / self.theta[0]*np.sqrt(self.theta[1]), X2 / self.theta[0]*np.sqrt(self.theta[1]), metric='sqeuclidean')
            K = (1 + dists)**(-self.theta[1])
        return K
    

class Periodic():
    def __init__(self, frequency=1.0, frequency_bounds=(1e-5, 1e5)):
        self.theta = [frequency]
        self.bounds = [frequency_bounds]
    def __call__(self, X1, X2=None):
        if X2 is None:
            # K(X1, X1) is symmetric so avoid redundant computation using pdist.
            dists = pdist(X1, lambda xi, xj: np.dot(np.sin(self.theta[0] * np.pi * (xi - xj)).T, 
                np.sin(self.theta[0] * np.pi * (xi - xj))))
            K = np.exp(-dists)
            K = squareform(K)
            np.fill_diagonal(K, 1)
        else:
            dists = cdist(X1, X2, lambda xi, xj: np.dot(np.sin(self.theta[0] * np.pi * (xi - xj)).T, 
                np.sin(self.theta[0] * np.pi * (xi - xj))))
            K = np.exp(-dists)
        return K