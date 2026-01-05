import unittest
import numpy as np
import os
import sys

# Add parent dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from lcfv import compute_dense_sift, gmm_em, fv_train, fv_encode, solve_lcfv1, solve_lcfv2

class TestSIFT(unittest.TestCase):
    def test_output_shape(self):
        img = np.random.rand(64, 64)
        descs = compute_dense_sift(img, step_size=4, bin_size=4)
        self.assertEqual(descs.shape[0], 128)
        self.assertGreater(descs.shape[1], 0)
        
    def test_rgb(self):
        img = np.random.rand(32, 32, 3)
        descs = compute_dense_sift(img)
        self.assertEqual(descs.shape[0], 128)
        
    def test_zero(self):
        img = np.zeros((32, 32))
        descs = compute_dense_sift(img)
        self.assertTrue(np.all(descs == 0), "Zero image should produce zero descriptors")

class TestFisher(unittest.TestCase):
    def test_gmm(self):
        np.random.seed(42)
        X = np.random.randn(2, 300)
        K = 3
        w, mu, sigma = gmm_em(X, K, max_iters=5)
        self.assertEqual(len(w), K)
        self.assertEqual(mu.shape, (2, K))
        self.assertEqual(sigma.shape, (2, K))
        self.assertFalse(np.any(np.isnan(mu)))
        
    def test_fv_encode_shape(self):
        D = 10
        K = 4
        X = np.random.rand(D, 50)
        w = np.ones(K)/K
        mu = np.random.rand(D, K)
        sigma = np.ones((D, K))
        
        fv = fv_encode(X, w, mu, sigma)
        self.assertEqual(fv.shape, (2*D*K, 1))
        
    def test_pca(self):
        X = np.random.randn(10, 100)
        _, _, _, pca_transform, _ = fv_train(X, 1, pca_dim=5)
        self.assertEqual(pca_transform.shape, (5, 10))

class TestLCFV(unittest.TestCase):
    def test_lcfv(self):
        dim = 20
        N = 50
        G = np.random.randn(dim, N)
        labels = np.random.randint(0, 2, N)
        
        C = (labels[:, None] == labels[None, :]).astype(float)
        alpha = 10
        
        # Test LCFV1
        M1, W1 = solve_lcfv1(G, C, alpha)
        # M1 acts on G
        feat1 = np.dot(M1, G)
        self.assertEqual(feat1.shape, (dim, N))
        
        # Test LCFV2
        M2 = solve_lcfv2(G, C, alpha)
        feat2 = np.dot(M2, G)
        # size depends on rank of C which is random, but at least dim + rank
        self.assertEqual(feat2.shape[1], N)

if __name__ == '__main__':
    unittest.main()
