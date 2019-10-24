## @package teetool
#  This module contains the GaussianProcess class
#
#  See GaussianProcess class for more details

import teetool as tt
import numpy as np

from numpy.linalg import det, inv, pinv, cond

## GaussianProcess class evaluates an ensemble of trajectories as a Gaussian process
#
#  Such a Gaussian process has a mean and covariance, and expresses itself as an ellipse (2d) or ellipsoid (3d) at a constant variance
class GaussianProcess(object):

    ## The constructor of GaussianProcess
    #   @param self             object pointer
    #   @param cluster_data     trajectory data in specific format: a list of (x, Y), where x [npoints x 1] and Y [npoints x ndim]
    #   @param ngaus            number of Gaussians desired
    def __init__(self, cluster_data, ngaus):
        # Fit cluster_data in a [0, 1] domain
        outline = tt.helpers.get_cluster_data_outline(cluster_data)
        cluster_data_norm = tt.helpers.get_cluster_data_norm(cluster_data,
                                                             outline)

        ## normalised cluster_data
        self._cluster_data = cluster_data_norm
        ## original outline
        self._outline = outline
        ## number of Gaussians after modelling
        self._ngaus = ngaus
        ## dimensionality of trajectory data
        self._ndim = tt.helpers.getDimension(cluster_data)

    ## obtain vectors to multiply normalised values with, allows for a transformation back to the actual values from the normalised ones.
    #  @param self  The object pointer.
    #  @return      M, vector with minimum values [ngaus*ndim x 1]
    #  @return      D, vector with difference values [ngaus*ndim x 1]
    def _outline2vectors(self):
        M_list = [] # list with minimum values [M]
        D_list = [] # list with difference [D]
        for d in range(self._ndim):
            xmin = self._outline[d*2+0]
            xmax = self._outline[d*2+1]

            m1 = np.ones(shape=(self._ngaus, 1)) * (xmin)
            M_list.append(m1)

            d1 = np.ones(shape=(self._ngaus, 1)) * (xmax - xmin)
            D_list.append(d1)

        M = np.concatenate(M_list, axis=0) # vector
        D = np.concatenate(D_list, axis=0) # vector

        return (M, D)

    ## returns the mu_y, sig_y vector to the original dimensions using the outline
    #  @param self  The object pointer.
    #  @param mu_y   mean vector [ngaus*ndim x 1]
    #  @param sig_y  covariance matrix [ngaus*ndim x ngaus*ndim]
    def _norm2real(self, mu_y, sig_y):
        (M, D) = self._outline2vectors()

        D_diag = np.diagflat(D ** 2)

        mu_y_real = mu_y * D + M
        sig_y_real = sig_y @ D_diag

        return mu_y_real, sig_y_real

    ## models the trajectory data via re-sampling, ignoring noise, missing data, trends, etc. Quick method only suitable for high-quality data
    #  @param self  The object pointer.
    #  @return mu_y mean vector [ngaus*ndim x 1]
    #  @return sig_y covariance matrix [ngaus*ndim x ngaus*ndim]
    #  @return cc mean [ndim x 1] in ngaus cells
    #  @return cA covariance [ndim x ndim] in ngaus cells
    def model_by_resampling(self):
        # extract
        cluster_data = self._cluster_data
        ngaus = self._ngaus
        mdim = self._ndim

        # predict these values
        xp = np.linspace(0, 1, ngaus)

        yc = []  # list to put trajectories

        for (xn, Yn) in cluster_data:
            # array to fill
            yp = np.empty(shape=(ngaus, mdim))

            for d in range(mdim):
                ynd = Yn[:, d]
                yp[:, d] = np.interp(xp, xn, ynd)

            # single column
            yp1 = np.reshape(yp, (-1, 1), order='F')

            yc.append(yp1)

        # compute values
        ntraj = len(yc)  # number of trajectories

        # obtain average [mu]
        mu_y = np.zeros(shape=(mdim*ngaus, 1))

        for yn in yc:
            mu_y += yn

        mu_y = (mu_y / ntraj)

        # obtain standard deviation [sig]
        sig_y_sum = np.zeros(shape=(mdim*ngaus, mdim*ngaus))

        for yn in yc:
            sig_y_sum += (yn - mu_y) @ (yn - mu_y).T

        sig_y = sig_y_sum / ntraj

        # convert to original values
        mu_y, sig_y = self._norm2real(mu_y, sig_y)

        # convert to cells
        (cc, cA) = self._getGMMCells(mu_y, sig_y, self._ngaus)

        return (mu_y, sig_y, cc, cA)

    ## models the trajectory data via maximum likelihood. It uses the basis function as specified to handle missing data, however, noise per trajectory has no influence on the parameter estimation. A suitable method in the absence of noise and known shape of trajectories.
    #  @param self  The object pointer.
    #  @param type_basis see Basis class for input
    #  @param nbasis see Basis class for input
    #  @return mu_y mean vector [ngaus*ndim x 1]
    #  @return sig_y covariance matrix [ngaus*ndim x ngaus*ndim]
    #  @return cc mean [ndim x 1] in ngaus cells
    #  @return cA covariance [ndim x ndim] in ngaus cells
    def model_by_ml(self, type_basis, nbasis):
        # extract
        cluster_data = self._cluster_data
        ngaus = self._ngaus

        ndim = self._ndim
        ntraj = len(cluster_data)

        # create a basis
        basis = tt.basis.Basis(type_basis, nbasis, ndim)

        wc = []

        for i, (xn, Y) in enumerate(cluster_data):
            yn = np.reshape(Y, newshape=(-1,1), order='F')
            Hn = basis.get(xn)
            wn = pinv(Hn) @ yn
            wn = wn
            wc.append(wn)

        # obtain average [mu]
        mu_w = np.zeros(shape=(ndim*nbasis, 1))

        for wn in wc:
            mu_w += wn

        mu_w = mu_w / ntraj

        # obtain standard deviation [sig]
        sig_w_sum = np.zeros(shape=(ndim*nbasis, ndim*nbasis))

        for wn in wc:
            sig_w_sum += (wn - mu_w) @ (wn - mu_w).T

        sig_w = sig_w_sum / ntraj

        # predict these values
        xp = np.linspace(0, 1, ngaus)
        Hp = basis.get(xp)

        mu_y = Hp @ mu_w
        sig_y = Hp @ sig_w @ Hp.T

        # convert to original values
        mu_y, sig_y = self._norm2real(mu_y, sig_y)

        (cc, cA) = self._getGMMCells(mu_y, sig_y, self._ngaus)

        return (mu_y, sig_y, cc, cA)

    ## models the trajectory data via expectation maximization. It uses the basis function as specified to handle missing data, and, when noisy data is detected within a trajectory, the global trend, as learned, takes over. A suitable method in the presence of noise or an unknown shape of trajectories -- the latter as different models can be compared via likelihood
    #  @param self  The object pointer.
    #  @param type_basis see Basis class for input
    #  @param nbasis see Basis class for input
    #  @param maximum_iterations maximum allowed number of evaluations till convergence
    #  @return mu_y mean vector [ngaus*ndim x 1]
    #  @return sig_y covariance matrix [ngaus*ndim x ngaus*ndim]
    #  @return cc mean [ndim x 1] in ngaus cells
    #  @return cA covariance [ndim x ndim] in ngaus cells
    def model_by_em(self, type_basis, nbasis, maximum_iterations=2001):
        # extract
        cluster_data = self._cluster_data
        ngaus = self._ngaus

        ndim = self._ndim
        ntraj = len(cluster_data)

        Mstar = 0
        for (xn, Yn) in cluster_data:
            Mstar += np.size(xn)

        # create a basis
        basis = tt.basis.Basis(type_basis, nbasis, ndim)

        # from cluster_data to cell structure
        yc, Hc = self._from_clusterdata2cells(cluster_data, basis)

        # hardcoded parameters
        MAX_ITERATIONS = maximum_iterations  # maximum number of iterations
        CONV_LIKELIHOOD = 1e-3  # stop convergence
        # min_eig = 10**-6  # minimum eigenvalue (numerical trick)

        # initial variables
        BETA_EM = 1000.
        mu_w = np.zeros(shape=(nbasis*ndim, 1))
        sig_w = np.eye(nbasis*ndim)
        sig_w_inv = inv(sig_w)

        loglikelihood_previous = np.inf

        for i_iter in range(MAX_ITERATIONS):

            # Expectation (54) (55)
            (Ewc, Ewwc) = self._Ewc_Ewwc(yc, Hc, mu_w, sig_w_inv, BETA_EM)

            #  Maximization :: (56), (57)

            # E [ MU ]
            mu_w = self._E_mu(Ewc)

            # E [ SIGMA ]
            sig_w = self._E_sigma(mu_w, yc, Hc, Ewc, Ewwc)

            # pre-calculate inverse
            sig_w_inv = inv(sig_w)

            # E [BETA]
            BETA_EM = self._E_beta(yc, Hc, Ewc, Ewwc, ndim, Mstar)

            # ////  log likelihood ///////////

            # // ln( p(Y|w) - likelihood
            loglikelihood_pYw = self._L_pYw(yc,
                                            Hc,
                                            Ewc,
                                            Ewwc,
                                            ndim,
                                            Mstar,
                                            BETA_EM)

            # // ln( p(w) ) - prior
            loglikelihood_pw = self._L_pw(Ewc,
                                          Ewwc,
                                          mu_w,
                                          sig_w,
                                          sig_w_inv,
                                          ndim,
                                          nbasis)

            loglikelihood_pY = loglikelihood_pYw + loglikelihood_pw

            # // check convergence
            loglikelihood_diff = np.abs(loglikelihood_pY - loglikelihood_previous)

            if np.isfinite(loglikelihood_pY):
                # check
                if (loglikelihood_diff < CONV_LIKELIHOOD):
                    break
            else:
                # not a valid loglikelihood
                print("warning: not a finite loglikelihood")
                break

            # output
            #if (i_iter % 100 == 0):
            #    print("{0} {1} {2}".format(i_iter, loglikelihood_pY, min_eig))

            # store previous log_likelihood
            loglikelihood_previous = loglikelihood_pY

        # predict these values
        xp = np.linspace(0, 1, ngaus)
        Hp = basis.get(xp)

        mu_y = Hp @ mu_w
        sig_y = Hp @ sig_w @ Hp.T

        # convert to original values
        mu_y, sig_y = self._norm2real(mu_y, sig_y)

        (cc, cA) = self._getGMMCells(mu_y, sig_y, self._ngaus)

        return (mu_y, sig_y, cc, cA)

    def _from_clusterdata2cells(self, cluster_data, basis):
        """converts from cluster_data (xn, Yn) list, to cells

        Input:
            cluster_data -
            basis -

        Output:
            yc -
            Hc -
        """

        # prepare data
        yc = []
        Hc = []

        for (xn, Yn)  in cluster_data:
            # data
            yn = np.reshape(Yn, newshape=(-1,1), order='F')
            Hn = basis.get(xn)
            # add to list
            yc.append(yn)
            Hc.append(Hn)

        return (yc, Hc)


    def _Ewc_Ewwc(self, yc, Hc, mu_w, sig_w_inv, BETA_EM):
        """returns the expected values Ewc and Ewwc

        input:
            yc          - [points]
            Hc          - [Gram matrix]
            mu_w        - E[w]
            sig_w_inv   - 1 / E[ww]
            BETA_EM     - 1 / noise

        output:
            Ewc     - [E[wn]]
            Ewnwc   - [E[wnwn]]
        """

        ntraj = len(yc)

        Ewc = []
        Ewwc = []

        # Expectation (54) (55)
        for n  in range(ntraj):
            # data
            yn = yc[n]
            Hn = Hc[n]

            (Ewn, Ewnwn) = self._Ewn_Ewnwn(yn,
                                           Hn,
                                           mu_w,
                                           sig_w_inv,
                                           BETA_EM)

            # store
            Ewc.append(Ewn);
            Ewwc.append(Ewnwn);

        return (Ewc, Ewwc)


    def _Ewn_Ewnwn(self, yn, Hn, mu_w, sig_w_inv, BETA_EM):
        """returns the expected values Ewn and Ewnwn

        input:
            yn          - points
            Hn          - Gram matrix
            mu_w        - E[w]
            sig_w_inv   - 1 / E[ww]
            BETA_EM     - 1 / noise

        output:
            Ewn     - E[wn]
            Ewnwn   - E[wnwn]
        """

        # calculate S :: (50)
        Sn_inv = sig_w_inv + BETA_EM * (Hn.T @ Hn)
        Sn = inv(Sn_inv)

        Ewn = (Sn @ ( (BETA_EM * (Hn.T @ yn)) + (sig_w_inv @ mu_w) ) )

        # BISHOP (2.62)
        Ewnwn = Sn + Ewn @ Ewn.T

        return (Ewn, Ewnwn)

    def _E_mu(self, Ewc):
        """returns the expected value E [ MU ]

        Input:
            Ewc - list of expected values

        Output:
            mu_w - average of expected values
        """

        # total number of trajectories
        ntraj = len(Ewc)

        mu_w_sum = np.zeros_like(Ewc[0])

        for Ewn in Ewc:
            # sum
            mu_w_sum += Ewn

        mu_w = mu_w_sum / ntraj

        return mu_w

    def _E_sigma(self, mu_w, yc, Hc, Ewc, Ewwc):
        """return the expected variance E [ SIGMA ]
        this takes into account the measured data and the model

        Input:
            mu_w -
            yc -
            Hc -
            Ewc -
            Ewwc -

        Output:
            sig_w -

        """

        # total number of trajectories
        ntraj = len(yc)

        sig_w_sum = np.zeros_like(Ewwc[0])

        # E [ SIGMA ]
        # sig_w_sum = np.zeros((nbasis*ndim, nbasis*ndim));

        for n  in range(ntraj):
            # extract data
            yn = yc[n]
            Hn = Hc[n]
            Ewn = Ewc[n]
            Ewnwn = Ewwc[n]

            # sum
            SIGMA_n = Ewnwn - 2.*(mu_w @ Ewn.T) + mu_w @ mu_w.T
            sig_w_sum += SIGMA_n

        sig_w = sig_w_sum / ntraj

        return sig_w

    def _E_beta(self, yc, Hc, Ewc, Ewwc, ndim, Mstar):
        """returns the expected noise parameter"""

        ntraj = len(yc)

        # E [BETA]
        BETA_sum_inv = 0.;

        for n  in range(ntraj):
            # extract data
            yn = yc[n]
            Hn = Hc[n]
            Ewn = Ewc[n]
            Ewnwn = Ewwc[n]

            BETA_sum_inv += yn.T @ yn - 2.*(yn.T @ (Hn @ Ewn)) + np.trace((Hn.T @ Hn) @ Ewnwn)

        BETA_EM = (ndim*Mstar) / BETA_sum_inv

        return BETA_EM

    def _L_pYw(self, yc, Hc, Ewc, Ewwc, ndim, Mstar, BETA_EM):
        """returns ln( p (Y|w) )
        likelihood of data, given the parameters"""

        ntraj = len(yc)

        loglikelihood_pYw_sum = 0.;

        for n  in range(ntraj):
            # extract data
            yn = yc[n]
            Hn = Hc[n]
            Ewn = Ewc[n]
            Ewnwn = Ewwc[n]

            # loglikelihood_pYw_sum = loglikelihood_pYw_sum + ((yn.')*yn - 2*(yn.')*(Hn*Ewn) + trace(((Hn.')*Hn)*Ewnwn));
            loglikelihood_pYw_sum += yn.T @ yn - 2.*(yn.T @ (Hn @ Ewn)) + np.trace((Hn.T @ Hn) @ Ewnwn)

        #  loglikelihood_pYw =  + ((Mstar*D) / 2) * log(2*pi) - ((Mstar*D) / 2) * log( BETA_EM ) + (BETA_EM/2) * loglikelihood_pYw_sum;
        loglikelihood_pYw = (Mstar*ndim / 2.) * np.log(2.*np.pi) - (Mstar*ndim / 2.) * np.log(BETA_EM) + (BETA_EM / 2.) * loglikelihood_pYw_sum

        return loglikelihood_pYw

    def _L_pw(self, Ewc, Ewwc, mu_w, sig_w, sig_w_inv, ndim, nbasis):
        """returns ln( p(w) )
        likelihood of parameters, before seeing the data"""

        # test conditioning sig_w
        if cond(sig_w) > 0:
            return 1e4

        ntraj = len(Ewc)

        loglikelihood_pw_sum = 0.;

        for n  in range(ntraj):
            # extract data
            Ewn = Ewc[n]
            Ewnwn = Ewwc[n]

            # loglikelihood_pw_sum = loglikelihood_pw_sum + trace( (LAMBDA_EM)*( Ewnwn - 2*MU_EM*(Ewn.') + (MU_EM*(MU_EM.')) ) );
            loglikelihood_pw_sum += np.trace(sig_w_inv*(Ewnwn - 2.*mu_w @ Ewn.T + mu_w @ mu_w.T))

        # loglikelihood_pw = + ((N*J*D) / 2) * log(2*pi) + (N/2) * ln_det_Sigma + (1/2) * loglikelihood_pw_sum;
        loglikelihood_pw = (ntraj*nbasis*ndim/2.)*np.log(2*np.pi) + (ntraj/2.)*np.log(det(sig_w)) + (1./2.)*loglikelihood_pw_sum

        return loglikelihood_pw

    def _getGMMCells(self, mu_y, sig_y, ngaus):
        """
        return Gaussian Mixture Model (GMM) in cells
        """

        cc = []
        cA = []

        for m in range(ngaus):
            # single cell
            (c, A) = self._getMuSigma(mu_y, sig_y, m, ngaus)

            # check for singularity
            A = tt.helpers.nearest_spd(A)

            cc.append(c)
            cA.append(A)

        return (cc, cA)


    def _getMuSigma(self, mu_y, sig_y, npoint, ngaus):
        """
        returns (mu, sigma)
        """
        # mu_y [DM x 1]
        # sig_y [DM x DM]
        D = self._ndim

        # check range
        if ((npoint < 0) or (npoint >= ngaus)):
            raise ValueError("{0}, not in [0, {1}]".format(npoint, ngaus))

        c = np.empty(shape=(D, 1))
        A = np.empty(shape=(D, D))

        # select position
        for d_row in range(D):
            c[d_row, 0] = mu_y[(npoint+d_row*ngaus), 0]
            for d_col in range(D):
                A[d_row, d_col] = sig_y[(npoint+d_row*ngaus), (npoint+d_col*ngaus)]

        return (c, A)
