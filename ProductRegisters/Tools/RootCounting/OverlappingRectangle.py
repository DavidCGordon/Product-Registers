from itertools import groupby

#group tuples by first dimension, and add an empty "0" layer
def preprocess(rectangle_list):
    output = sorted(rectangle_list, key = lambda x: x[0], reverse=True)
    output = [(k, [x[1:] for x in g]) for k, g in groupby(output, key = (lambda x: x[0]))]
    output.append((0,[]))
    return output

# check if tuple 1 is <= than tuple 2 in every index
# this indicates rectangle A is contained entirely in rectangle B
def rect_sset(rect_A,rect_B):
    return all(a <= b for a,b in zip(rect_A,rect_B))

#insert tuples from a into b, if they are still needed.
def insert_tuples(upper_layer, lower_layer):
    still_needed = []
    for upper_rectangle in upper_layer:
        needed = True

        # if any rectangle in the lower layer needs includes this one:
        # it is not needed for further area calculations
        for lower_rectangle in lower_layer:
            if rect_sset(upper_rectangle,lower_rectangle):
                needed = False
                break
        if needed:
            still_needed.append(upper_rectangle)

    #add the rectangles which are still needed:
    lower_layer += still_needed

#recursive call to solve (grows poorly with dim, but MUCH better with
#number of rectangles: significantly faster/more useable.
def _solve_rec(dimension,rectangle_list):
    # base case:
    if dimension == 1:
        return max(x[0] for x in rectangle_list)
    
    total_area = 0
    processed_list = preprocess(rectangle_list)

    for layer in range(len(processed_list)-1):

        # use the current dimension to get the width of the layer
        layer_width = processed_list[layer][0] - processed_list[layer+1][0]
        # solve the subproblem 1 dimension down to get the cross-sectional area
        cross_section_area = _solve_rec(dimension-1, processed_list[layer][1])
        # multiply to get the volume of the cross-section and add to total area
        total_area += layer_width * cross_section_area

        # insert any still needed rectangles into the next layer:
        # this ensures the cross-sectional area remains correct.
        insert_tuples(processed_list[layer][1], processed_list[layer+1][1])

    return total_area

#wrapper to preprocess dimension for ease of use.
def rectangle_solve(rectangle_list):
    dimension = len(rectangle_list[0])
    return _solve_rec(dimension,rectangle_list)