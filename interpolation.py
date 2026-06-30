

# =====================================================================
# Assignment: Polynomial Interpolation (Lagrange & Neville)
# Group Members:
#   1. [Student Name 1] - [ID 1]
#   2. [Student Name 2] - [ID 2]
#   3. [Student Name 3] - [ID 3]
# =====================================================================

def lagrange_interpolation(x_points, y_points, x_target):
    """
    Evaluates the interpolating polynomial at x_target using Lagrange Interpolation.

    Parameters:
    x_points (list): List of X coordinates of the data points.
    y_points (list): List of Y coordinates of the data points.
    x_target (float): The point at which to interpolate.

    Returns:
    float: The interpolated value at x_target.
    """
    n = len(x_points)
    y_target = 0.0

    # Outer loop for summing all terms: P(x) = Sigma (y_i * L_i(x))
    for i in range(n):
        l_i = 1.0
        # Inner loop to calculate the basis polynomial L_i(x)
        for j in range(n):
            if i != j:
                # Lagrange formula product term
                l_i *= (x_target - x_points[j]) / (x_points[i] - x_points[j])

        # Add the current term to the total sum
        y_target += y_points[i] * l_i

    return y_target


def neville_interpolation(x_points, y_points, x_target):
    """
    Evaluates the interpolating polynomial at x_target using Neville's Algorithm.

    Parameters:
    x_points (list): List of X coordinates of the data points.
    y_points (list): List of Y coordinates of the data points.
    x_target (float): The point at which to interpolate.

    Returns:
    float: The interpolated value at x_target.
    """
    n = len(x_points)

    # Initialize the tableau (table) with zeros
    # Q[i][j] represents the polynomial value passing through points i to i+j
    Q = [[0.0] * n for _ in range(n)]

    # Base case: degree 0 polynomials are just the y values themselves
    for i in range(n):
        Q[i][0] = float(y_points[i])

    # Iteratively fill the Neville tableau column by column
    for j in range(1, n):  # j represents the degree of the polynomial
        for i in range(n - j):  # i represents the starting index of the points
            # Apply Neville's recursive formula
            numerator1 = (x_target - x_points[i + j]) * Q[i][j - 1]
            numerator2 = (x_target - x_points[i]) * Q[i + 1][j - 1]
            denominator = x_points[i] - x_points[i + j]

            Q[i][j] = (numerator1 - numerator2) / denominator

    # The final result is at the top right of the completed tableau
    return Q[0][n - 1]


def main():
    print("==================================================")

    # 1.1 Defining the table data points (Inputs)
    # Example data points from a known polynomial function (e.g., y = x^2)
    x_points = [1.0, 2.0, 3.0, 4.0]
    y_points = [1.0, 4.0, 9.0, 16.0]

    # 1.2 The target point where we want to find the estimated value (Input)
    x_target = 2.5

    # Display input data clearly
    print("Input Data Points configured in the system:")
    for i in range(len(x_points)):
        print(f"  Point {i + 1}: ({x_points[i]}, {y_points[i]})")
    print(f"Target point to estimate (X): {x_target}")
    print("==================================================")

    # Executing both interpolation methods
    result_lagrange = lagrange_interpolation(x_points, y_points, x_target)
    result_neville = neville_interpolation(x_points, y_points, x_target)

    # 1.3 Displaying the final results with appropriate descriptive messages
    print("Execution Results:")
    print(f"  -> Lagrange Interpolation Result : {result_lagrange}")
    print(f"  -> Neville's Algorithm Result     : {result_neville}")
    print("--------------------------------------------------")
    print("Note: Since the points are from y = x^2, the exact")
    print("theoretical value for X = 2.5 should be exactly 6.25.")
    print("==================================================")


if __name__ == "__main__":
    main()