""" illuminutils.domain.observer_position - Observer position index from coordinates

Input:
    corrected_bbox : Hardcoded corrected bbox
    pixsiz : Hardcoded pixel size
    (observer_xy, observer_alt, observer_height) : Hardcoded observer position (xy, altitude and height)
    Hardcoded minimal altitude of the domain

Output:
    Observer's cell index (in text)
    Vertical resolution at observer height (in texte)

"""
import math

if __name__ == '__main__':
    
    # Hardcoded inputs begin here
    corrected_bbox = [239981.58663178002,4324182.311412051,639981.58663178,4624182.311412051]
    pixsiz = 1000
    observer_xy = (239981.58663178002+200000,4324182.311412051+150000)
    observer_alt = 659.0 # From DEM
    observer_height = 1 # Observer height above the ground
    minimal_alt = 211
    # End of hardcoded inputs
    
    # Cell height are given in Illumina code
    cell_height = [0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,
            2.68,3.24,3.92,4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,
            25.94,31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
            171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,778.66,
            940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,3530.84,
            4265.26,5152.44]

    observer_cell_x = math.floor((observer_xy[0] - corrected_bbox[0])/pixsiz)
    observer_cell_y = math.floor((observer_xy[1] - corrected_bbox[1])/pixsiz)

    observer_delta_z = observer_alt + observer_height - minimal_alt

    # To find position of observer in z, we use the following loop:
    # First, we set base cumulative height to 0 meter and
    # observer z position to 0.
    # For each iteration, we check if cumulative height is now higher than
    # observer relative height. If so, we have found z position.
    # else,  we increase z index and we return to the beginning of loop
    cumulative_height = 0
    observer_cell_z = 0
    for case in range(len(cell_height)):
        cumulative_height += cell_height[case]
        if observer_delta_z < cumulative_height:
            break
        else:
            observer_cell_z = case

    print """Observer position in cell is (%i, %i, %i)""" % \
            (observer_cell_x, observer_cell_y, observer_cell_z)
    print """Vertical resolution at observer height is %.3f""" % \
            (cell_height[observer_cell_z],)
