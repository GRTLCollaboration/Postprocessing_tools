
# copy of GRChombo's C++ Integration Methods class

# A class to store and return the weights associated to a Newton-Cotes formula
# for numerical integration/quadrature which can be closed (i.e. includes the
# endpoints) or open (does not include the end points). This is used by
# SurfaceExtraction for integration over extraction surfaces.
class IntegrationMethod:
    ### public
    def __init__(self, weights, is_closed = True):
        self.weights = weights
        self.num_weights = len(weights)
        self.is_closed = is_closed
        assert(self.num_weights > 0);


    # Checks that this integration method is suitable given the number of points and periodicity
    def isValid(self, num_points, is_periodic):
        if self.is_closed and not is_periodic:
            return (num_points % self.num_weights == 1 or self.num_weights == 1);
        else:
            return (num_points % self.num_weights == 0);

    # Returns the weight for a point with given index
    def weight(self, index, total_points, is_periodic):
        weight_index = index % self.num_weights;
        endpoint = (index == 0 or index == total_points - 1) and not is_periodic;
        # if this is a closed formula, not a geometry endpoint but at the edge
        # of the formula, need to double the weight as this is how Newton-Cotes
        # formulae are combined.
        if self.is_closed and not endpoint and weight_index == 0:
            return 2. * self.weights[weight_index];
        else: # otherwise we just use the weight from the formula
            return self.weights[weight_index];

# define some static IntegrationMethods here
# Closed methods:
trapezium = IntegrationMethod([0.5]);
simpson = IntegrationMethod([0.3333333333333333, 1.3333333333333333]);
simpson38 = IntegrationMethod([0.375, 1.125, 1.125]);
boole = IntegrationMethod([0.3111111111111111, 1.4222222222222222,
                           0.53333333333333, 1.4222222222222222]);

# Open Methods:
midpoint = IntegrationMethod([1.0], False);
milne_regularized = IntegrationMethod([1.125, 0.75, 1.125], False);
open_3rd_order = IntegrationMethod([1.0833333333333333, 0.9166666666666666,
                                    0.9166666666666666, 1.0833333333333333],
                                   False);
open_4th_order = IntegrationMethod(
    [1.1935763888888888, 0.4340277777777778, 1.7447916666666667,
     0.4340277777777778, 1.1935763888888888],
    False);
