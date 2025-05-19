import numpy as np


class Controller:
    def __init__(self, kp, ki, td, kaw, output_limit):


        """
        PI(D) controller with anti-windup

        Key idea is to add lots of different functionalities to the controller to be more easely utilized in simulation.


        """


        self.kp = kp  # Proportional gain
        # self.ti = ti  # Integral time constant
        self.ki = ki
        self.td = td  # Derivative time constant
        self.kaw = kaw  # Anti-windup gain
        self.output_limit = output_limit  # Max controller output (e.g., DC-link limit)

        self.integral = 0.0
        self.last_output = 0.0
        self.prev_error = 0.0

        self.P_term = 0.0
        self.I_term = 0.0
        self.D_term = 0.0

    def compute(self, error, dt, time):
        # PID control
        self.integral += error * dt

        self.P_term = self.kp * error
        self.I_term = (self.ki) * self.integral
        self.D_term = (self.kp * self.td) * (error - self.prev_error) / dt

        u_desired = self.P_term + self.I_term + self.D_term

        u_out = np.clip(u_desired, -self.output_limit, self.output_limit)

        # Anti-windup correction
        aw_term = self.kaw * (u_desired - u_out)
        self.integral += aw_term * dt

        # Store output and prev error for next iteration
        self.last_output = u_out
        self.prev_error = error

        return u_out

# endregion
