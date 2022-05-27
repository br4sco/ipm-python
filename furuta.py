import math
import traceback
import logging
import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


class FurutaODE:
    """A class representing the Furuta pendulum

    Static Attributes
    ----------
    integration_method : str
        The integration method used when simulating the model

    ys : str list
        The output variable names where:
        t [s] is the time.

        theta [rad] is the angle between z-axis and the pendulum. A theta = 0
        indicates that the pendulum is standing up vertically, and theta = pi
        that the pendulum is hanging vertically.

        phi [rad] is the angle between the arm and x-axis.

        dthetadt [rad/s] and dphidt [rad/s] are the first derivatives of theta
        and phi with respect to t, respectivly. Angles are given in radians.

    ys_to_idx : { str : int }
        A dictionary mapping output variables names to output variable indexed
        in the state vector q. The variable t has index -1 and is not part of q.

    Attributes
    ----------
    q : numpy.array([float])
        The state vector of dependent variables

    t : float
        The free variable time

    wrap_angles : Bool
        If angles are confined to the interval [0, 2pi]

    Methods
    -------

    IPM interface methods
    ---------------------
    init(...)
        Initializes the model and simulation

    trans(...)
        Transitions the simulation a finite time step

    output(...)
        Outputs the current values of the output variables

    Helper methods
    ---------------
    plot(...)
        Computes a dense solution to the model over a time-interval and then
        plots the results
    """

    integration_method = "BDF"
    ys = ["theta", "phi", "dthetadt", "dphidt", "t"]
    ys_to_idx = {"theta": 0, "phi": 1, "dthetadt": 2, "dphidt": 3, "t": -1}

    def _check_ys(ys):
        if not all(y in FurutaODE.ys for y in ys):
            raise ValueError(f"ys expected to be a subset of {FurutaODE.ys}")

    def f(self, t, q, u, a, b, c, d2, K, J, b1, b2):
        """Implementation based on equation (2.10) in
        'Furuta’s Pendulum: A Conservative Nonlinear Model for Theory and
         Practise', J. Á. Acosta.
        """
        q1 = q[0]
        q2 = q[1]
        dq1 = q[2]
        dq2 = q[3]
        M = J * np.array(
            [
                [1.0, a * np.cos(q1)],
                [a * np.cos(q1), b + np.sin(q1) ** 2.0],
            ]
        )
        C = J * np.array(
            [
                [b1, -dq2 * np.sin(q1) * np.cos(q1)],
                [
                    -a * dq1 * np.sin(q1) + dq2 * np.sin(q1) * np.cos(q1),
                    dq1 * np.sin(q1) * np.cos(q1) + b2,
                ],
            ]
        )
        dvdq = J * np.array([-(d2) * np.sin(q1), 0.0])
        dqdt = np.array([dq1, dq2, 0.0, 0.0])
        dqdt[2:] += np.linalg.solve(
            M, np.array([0.0, K * u]) - C @ q[2:] - dvdq
        )
        return dqdt

    # def f_hamiltonian(self, t, q, u, a, b, c, d2, K, J):
    #     """Implements equaiton (2.22) in
    #     'Furuta’s Pendulum: A Conservative Nonlinear Model for Theory and
    #      Practise', J. Á. Acosta.
    #     """
    #     q1 = q[0]
    #     q2 = q[1]
    #     p1 = q[2]  # note that this is not dthetadt
    #     p2 = q[2]  # note that this it not dphidt
    #     det = J * (b + (np.sin(q1) ** 2) * p1 - (a * np.cos(q1)) ** 2)
    #     return np.array(
    #         [
    #             ((b + np.sin(q1) ** 2) * p1 - a * p2 * np.cos(q1)) / det,
    #             (p2 - a * p1 * np.cos(q1)) / det,
    #             (
    #                 J
    #                 * (
    #                     (b * np.sin(q1) ** 2) * p1**2
    #                     + p2**2
    #                     - 2 * a * p1 * p2 * np.cos(q1)
    #                 )
    #                 * ((1 + a**2) * np.sin(2 * q1))
    #             )
    #             / (2 * det**2)
    #             + J * d2 * np.sin(q1)
    #             - (p1**2 * np.sin(2 * q1) + 2 * a * p1 * p2 * np.sin(q1))
    #             / (2 * det),
    #             J * c * d2 * u,
    #         ]
    #     )

    def __init__(self, wrap_angles=True):

        self.q = None
        self.t = None
        self.wrap_angles = wrap_angles
        self._args = None

    def _solve(self, time, u, rtol, atol):
        tend = self.t + time
        return solve_ivp(
            fun=self.f,
            t_span=(self.t, tend),
            y0=self.q,
            args=(u, *self._args),
            rtol=rtol,
            atol=atol,
            dense_output=True,
            method=FurutaODE.integration_method,
        )

    def _check_init(self):
        if self.t is None or self.q is None or self._args is None:
            raise Exception("Not initialized, call init(...)")

    def init(
        self,
        t0=0.0,
        theta0=math.pi,
        phi0=0.0,
        dthetadt0=0.0,
        dphidt0=0.0,
        m=1.0,
        l=1.0,
        r=1.0,
        g=1.0,
        J=1.0,
        Ja=1.0,
        K=1.0,
        b1=0.01,
        b2=0.01,
    ):
        """Initializes the model and simulation

        Parameters
        ----------
        t0 : float, optional
            The initial time [s] (default is 0.0)

        theta0 : float, optional
            The inital value of theta [rad] (default is π)

        phi0 : float, optional
            The inital value of phi [rad] (default is 0.0)

        dthetadt0 : float, optional
            The inital value of dthetadt [rad] (default is 0.0)

        dphidt0 : float, optional
            The inital value of dtphidt [rad] (default is 0.0)

        m : float, optional
            The mass [kg] of the the pendulum (default is 1.0)

        l : float, optional
            Half the pendulum arm length [m] (default is 1.0)

        r : float, optional
            Radius of the arm [m] (default is 1.0)

        g : float, optional
            Acceleration of gravity [m / s^2] (default is 1.0)

        J : float, optional
            Moment of inertia of the pendulum w.r.t. the pivot [kg m^2]
            (default is 1.0)

        Ja : float, optional
            Moment of inertia of the arm and the motor [kg m^2] (default is
            1.0)

        K : float, optional
            Constant torque [N m / v] (default is 1.0)

        b1 : float, option
            Viscous damping coefficient in motor joint (default is 0.01)

        b2 : float, option
            Viscous damping coefficient in pivot joint (default is 0.01)
        """

        a = m * r * l / J
        b = (Ja + m * r**2.0) / J
        c = K / (m * g * l)
        d2 = m * g * l / J

        self._args = (a, b, c, d2, J, K, b1, b2)
        self.t = t0

        self.q = np.array(
            [
                theta0,
                phi0,
                dthetadt0,
                dphidt0,
            ]
        )
        if self.wrap_angles:
            self.q[:2] %= 2 * math.pi

    def trans(self, u, step, rtol=1e-4, atol=1e-6):
        """Transitions the simulation a finite time step

        Parameters
        ----------
        u : float
            Control input (controls the torque output of the motor)

        step : float
            The length [s] of this transition

        rtol : float, optional
            The relative tolarance of the underlying ode solver (default is 1e-4)

        atol : float, optional
            The absolute tolarance of the underlying ode solver (default is 1e-6)

        Raises
        ------
        Exception
            If trans(...) is called before init(...)

        Returns
        -------
        Boolean
            If the transition is successful. If the transition is not successful
            the internal state is not updated.
        """

        self._check_init()
        tend = self.t + step

        try:
            sol = solve_ivp(
                fun=self.f,
                t_span=(self.t, tend),
                y0=self.q,
                t_eval=[tend],
                args=(u, *self._args),
                rtol=rtol,
                atol=atol,
                method=FurutaODE.integration_method,
            )

            if sol.success:
                self.q = sol.y.flatten()
                if self.wrap_angles:
                    self.q[:2] %= 2 * math.pi
                self.t = sol.t[0]

            return sol.success

        except Exception as err:
            logging.error(
                f"An exception occurred during solving for input {u}: {err}"
            )
            return False

    def output(self, ys=["theta"]):
        """Outputs the current values of the output variables

        Parameters
        ----------
        ys : list str, optional
            The output variable names to get the output from. See documentation
            for the attribute ys (defaults to [\"theta\"])

        Raises
        ------
        ValueError
            If a requested output variable name does not exist

        Returns
        dict str float
            The output values for the requested output variables
        """

        FurutaODE._check_ys(ys)

        out = dict()
        for y in ys:
            if y == "t":
                out[y] = self.t
            else:
                out[y] = self.q[FurutaODE.ys_to_idx[y]]

        return out

    def plot(self, time=2.0, u=0.0, rtol=1e-4, atol=1e-6, ys=["theta"]):
        """Computes a dense solution to the model over a time-interval and then
        plots the results

        Parameters
        ---------
        time : float, optional
            The length of the simulation interval (default is 2.0)

        samples : float, optional
            Control input (controls the torque output of the motor,
            default is 0.0)

        rtol : float, optional
            The relative tolarance of the underlying ode solver (default is 1e-4)

        atol : float, optional
            The absolute tolarance of the underlying ode solver (default is 1e-6)

        ys : list str, optional
            The output variable names to plot the output from. See documentation
            for the attribute ys (defaults to [\"theta\"])

        Raises
        ------
        ValueError
            If a requested output variable name does not exist

        Exception
            If plot(...) is called before init(...)
        """

        FurutaODE._check_ys(ys)
        self._check_init()

        idx = [FurutaODE.ys_to_idx[y] for y in ys]
        if -1 in idx:
            idx.remove(-1)

        sol = self._solve(time, u, rtol, atol)
        h = 0.01
        n = int((time) / h)
        t = np.linspace(self.t, self.t + time, n)

        z = sol.sol(t)
        if self.wrap_angles:
            z[:2, :] %= 2 * math.pi

        legend = np.array(
            [r"$\theta$", r"$\phi$", r"$\dot{\theta}$", r"$\dot{\phi}$"]
        )
        plt.plot(t, z[idx, :].T)
        plt.xlabel("t")
        plt.legend(
            legend[idx],
            shadow=True,
        )
        plt.title("Furuta Pendulum")
        plt.show()


##################################
# Functional style IPM interface #
##################################


def init(
    model,
    t0=0.0,
    theta0=math.pi,
    phi0=0.0,
    dthetadt0=0.0,
    dphidt0=0.0,
    m=1.0,
    l=1.0,
    r=1.0,
    g=1.0,
    J=1.0,
    Ja=1.0,
    K=1.0,
):
    """See documention for FurutaODE.init(...).

    An instantiation of FurutaODE serves as the model parameter.

    """

    model = copy.deepcopy(model)
    model.init(t0, theta0, dthetadt0, dphidt0, m, l, r, g, J, Ja, K)
    return model


def trans(state, u, step, rtol=1e-4, atol=1e-6):
    """See documention for FurutaODE.trans(...).

    An instantiation of FurutaODE serves as the state parameter.

    """

    new_state = copy.deepcopy(state)
    if new_state.trans(u, step, rtol, atol):
        return new_state
    return state


def output(state, ys=["theta"]):
    """See documention for FurutaODE.output(...).

    An instantiation of FurutaODE serves as the state parameter.

    """

    state = copy.deepcopy(state)
    return state.output(ys)
