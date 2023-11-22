import random
from scipy import optimize
from ROOT import TF2
import math

def fit_primary_vertex(tracks, beam, z_target):

    x0 = beam["beam.vx"]*z_target+beam["beam.px"]
    y0 = beam["beam.vy"]*z_target+beam["beam.py"]

    #define the function
    def distToPvSquared(parameters):
        xv, yv, zv = parameters
        dist2 = 0
        for track in tracks:
            dist2 += (track[0]*zv+track[2]-xv)**2+(track[1]*zv+track[3]-yv)**2
        dist2 += (beam["beam.vx"]*zv+beam["beam.px"]-xv)**2+(beam["beam.vy"]*zv+beam["beam.py"]-yv)**2
        return dist2

    # Initial guess
    initial_guess = [x0, y0 , z_target]
    # Optimize the function using the BFGS method
    result = optimize.minimize(distToPvSquared, initial_guess, method='BFGS', tol=1e-5)
    
    # Retrieve the results
    return result.x

def define_PV_selection(residuals, conf_file, type = "Pb"):

    tf2 = TF2("xygaus","xygaus",-0.20,0.20,-0.2,0.2)
    residuals.Fit(tf2,"MR0")

    sigma_x = tf2.GetParameter(2)
    sigma_y = tf2.GetParameter(4)
    mean_x = tf2.GetParameter(1)
    mean_y = tf2.GetParameter(3)
    write_selection(conf_file, [sigma_x,sigma_y,mean_x,mean_y], type)


def write_selection(conf_file, parameters = [0,0,0,0], type = "Pb"):
    counter = 0
    with open(conf_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break

            if f"SIGMA_X_{type}" in line.strip():
                counter_SIGMA_X = counter
            if f"SIGMA_Y_{type}" in line.strip():
                counter_SIGMA_Y = counter
            if f"MEAN_X_{type}" in line.strip():
                counter_MEAN_X = counter
            if f"MEAN_Y_{type}" in line.strip():
                counter_MEAN_Y = counter
            counter += 1
    
    # Open the file in read mode and read all its lines into a list
    with open(conf_file, 'r') as file:
        lines = file.readlines()

    lines[counter_SIGMA_X] = f"SIGMA_X_{type} : {parameters[0]} \n"
    lines[counter_SIGMA_Y] = f"SIGMA_Y_{type} : {parameters[1]} \n"
    lines[counter_MEAN_X] = f"MEAN_X_{type} : {parameters[2]} \n"
    lines[counter_MEAN_Y] = f"MEAN_Y_{type} : {parameters[3]} \n"
    # Open the file in write mode and write the modified list to the file
    with open(conf_file, 'w') as file:
        file.writelines(lines)

        import math

def compute_dca(point, line_params):
    """
    Calculate the distance between a point and a line parameterized in 3D space.

    Parameters:
        point (tuple): A tuple representing the coordinates of the point (x, y, z).
        line_params (tuple): A tuple representing the parameters of the line (vx, vy, x0, y0).

    Returns:
        float: The distance between the point and the line in z, y, and z.
    """
    x1, y1, z1 = point
    vx, vy, x0, y0 = line_params

    numerator = (x1 - x0) * vx + (y1 - y0) * vy + z1
    denominator = vx**2 + vy**2 + 1

    t = numerator / denominator

    distance = [vx*t + x0 - x1, vy*t + y0 - y1, t - z1]
    return distance
