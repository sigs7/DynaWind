# pi_controller.py


# region PI Controller class
# class PIController_LAH:
#     def __init__(self, kp, ti):
#         self.kp = kp
#         self.ti = ti
#         self.integral = 0.0

#     def compute(self, error, dt):
#         self.integral += error * dt

#         # Anti-windup
#         if self.integral > 1.5:
#             self.integral = 1.5
#         if self.integral < -1.5:
#             self.integral = -1.5
        
#         return self.kp * error + (self.kp/self.ti) * self.integral


class PIController_LAH:
    def __init__(self, kp, ki, td, kaw, output_limit):
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


        # u_desired = self.kp*error + (self.kp / self.ti) * self.integral + (self.kp * self.td) * (error - self.prev_error) / dt
        u_desired = self.P_term + self.I_term + self.D_term

        # # Check if output is saturating
        # if u_desired > self.output_limit:
        #     u_actual = self.output_limit
        #     anti_windup_term = self.kaw * (u_actual - u_desired)
        #     self.integral -= anti_windup_term * dt  # Apply only during saturation
        # elif u_desired < -self.output_limit:
        #     u_actual = -self.output_limit
        #     anti_windup_term = self.kaw * (u_actual - u_desired)
        #     self.integral -= anti_windup_term * dt  # Apply only during saturation
        # else:
        #     u_actual = u_desired  # No saturation, no anti-windup applied

        # # Update internal variables
        # self.last_output = u_actual
        # self.prev_error = error

        # return u_actual

        return u_desired

# endregion
