# pi_controller.py


# region PI Controller class
class PIController_LAH:
    def __init__(self, kp, ti):
        self.kp = kp
        self.ti = ti
        self.integral = 0.0

    def compute(self, error, dt):
        self.integral += error * dt

        # Anti-windup
        if self.integral > 1.5:
            self.integral = 1.5
        if self.integral < -1.5:
            self.integral = -1.5
        
        return self.kp * error + (self.kp/self.ti) * self.integral

# endregion
