import matplotlib.pyplot as plt
import numpy as np
from numba import njit


@njit
def vector_from_angle(alpha_rad, beta_in):
    alpha_in = np.rad2deg(alpha_rad)
    lift = 0.0
    if 0.0 < alpha_in <= 20.0:
        lift = 1.0 - alpha_in / 20.0
    # if 160.0 < alpha_in <= 180.0:
    #     lift = (160.0 - alpha_in) / 20.0

    friction = 0.0
    if 20.0 < alpha_in <= 90.0 - beta_in:
        fraction_top = alpha_in - 10.0 - ((90.0 - beta_in) / 2.0)
        fraction_bottom = 45.0 - ((beta_in + 20.0) / 2.0)
        friction = 1.0 - abs(fraction_top / fraction_bottom)
    if 90.0 + beta_in < alpha_in < 160.0:
        fraction_top = alpha_in - 45.0 - ((beta_in + 160.0) / 2.0)
        fraction_bottom = 45.0 - ((beta_in + 20.0) / 2.0)
        friction = 1.0 - abs(fraction_top / fraction_bottom)

    return (lift, friction, 0.0)


@njit
def run():
    alpha = 0.01
    omega = 0.1
    epsilon = 0.0

    beta = 10.0

    m = 10.0 ** (-3.0)
    w = 0.03
    l = 0.3
    rho = 1.225
    g = 9.8
    C_d = 1.1

    dt = 0.000001
    I = (1.0 / 12.0) * m * (w ** 2.0)

    v_x = 0.1
    v_y = -0.1

    x = 0.0
    y = 0.0

    # Ensure num_steps is an integer for the array size
    num_steps = int(10.0 / dt)

    # Numba MUST use pre-allocated arrays, not standard Python lists
    x_s = np.zeros(num_steps)
    y_s = np.zeros(num_steps)
    alpha_s = np.zeros(num_steps)
    omega_s = np.zeros(num_steps)
    mode_vector_s = np.zeros((num_steps, 3))

    progress = 0.0

    i = 0
    while i < num_steps:
        epsilon = 0.0

        # Note: The string variable 'mode' was removed.
        # Numba cannot handle dynamic string assignments efficiently.

        if np.deg2rad(0.0) < alpha <= np.deg2rad(20.0):
            # Using tuples of floats instead of lists
            mode_vector = (1.0, 0.0, 0.0)
        elif np.deg2rad(20.0) < alpha <= np.deg2rad(90.0 - beta):
            mode_vector = (0.0, 1.0, 0.0)
        elif np.deg2rad(90.0 - beta) < alpha <= np.deg2rad(90.0 + beta):
            mode_vector = (0.0, 0.0, 1.0)
        elif np.deg2rad(90.0 + beta) < alpha <= 180.0:
            mode_vector = (0.0, 1.0, 0.0)
        else:
            mode_vector = (0.0, 0.0, 0.0)

        (lift_component, friction_component, _) = vector_from_angle(alpha, beta)
        distortion = 100000.0
        lift_component_modified_signless = 1.0 - (1.0 - abs(lift_component))**distortion
        if lift_component != 0.0:
            lift_component_modified_signed = lift_component_modified_signless*(lift_component/abs(lift_component))
        else:
            lift_component_modified_signed = 0.0
        mode_vector = (lift_component_modified_signed, 1.0 - (1.0 - friction_component)**distortion, 0.0)

        # Mode 1
        L = np.pi * rho * l * (v_x ** 2.0) * alpha
        M_LE = np.pi * alpha / 2.0
        torque = M_LE * (0.5 * rho * v_x ** 2.0) * (w * l) * w

        v_y += ((L / m) * dt) * mode_vector[0]
        epsilon += (torque / I) * mode_vector[0]

        # Mode 2
        v_perp_x = v_x * ((np.sin(alpha)) ** 2.0) - v_y * (np.sin(alpha)) * (np.cos(alpha))
        v_perp_y = v_y * ((np.cos(alpha)) ** 2.0) - v_x * (np.sin(alpha)) * (np.cos(alpha))

        # Replaced **(1/2) with **0.5
        v_perp = (v_perp_x ** 2.0 + v_perp_y ** 2.0) ** 0.5

        if v_perp != 0.0:
            v_perp_x_unit = v_perp_x / v_perp
            v_perp_y_unit = v_perp_y / v_perp
        else:
            v_perp_x_unit = 0.0
            v_perp_y_unit = 0.0

        k = 2.0

        v_x -= ((v_perp_x_unit * k * (v_perp ** 2.0)) * dt) * mode_vector[1]
        v_y -= ((v_perp_y_unit * k * (v_perp ** 2.0)) * dt) * mode_vector[1]

        if 0.01 * np.floor(100.0 * float(i) * dt) > progress:
            progress = 0.01 * np.floor(10.0 * float(i) * dt)
            print(int(100.0 * progress), "%, omega =", omega, ", alpha =", alpha)
            if omega < 0.0:
                print("OMEGA NEGATIVE")

        v_y -= g * dt
        draw_torque = ((1.0 / 64.0) * C_d * rho * (omega ** 2.0) * (w ** 4.0) * l)
        II = I
        epsilon_from_drag = -(((1.0 / 64.0) * C_d * rho * (omega * abs(omega)) * (w ** 4.0) * l) / I)
        epsilon += epsilon_from_drag

        change_in_omega = epsilon * dt
        omega += epsilon * dt

        x += v_x * dt
        y += v_y * dt
        alpha += omega * dt

        alpha = alpha % np.deg2rad(180.0)

        # Store in arrays using the index 'i'
        x_s[i] = x
        y_s[i] = y
        alpha_s[i] = alpha
        omega_s[i] = omega
        mode_vector_s[i, 0] = mode_vector[0]
        mode_vector_s[i, 1] = mode_vector[1]
        mode_vector_s[i, 2] = mode_vector[2]

        i += 1

    return (x_s, y_s, alpha_s, omega_s, mode_vector_s)


(x_s, y_s, alpha_s, omega_s, mode_vector_s) = run()

# Since mode_vector_s is a NumPy array, we use slicing instead of list comprehension
mode_vector_s_friction = mode_vector_s[:, 1]

print("Loop done")
print(alpha_s[1:100])
plt.scatter(x_s[1:-1:10], y_s[1:-1:10])

plt.show()