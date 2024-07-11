from sl_python.sl_2D import dep_search_1d

def test_dep_search_1d():
    xcr = 0
    i = 0
    fix1 = {
        "vx_e": 0.2,
        "vx_tmp": 0.2,
        "dth": 0.5,
        "dx": 0.02
    }
    
    fix2 = {
        "vx_e": -0.2,
        "vx_tmp": -0.2,
        "dth": 0.5,
        "dx": 0.02
    }
    
    fix3 = {
        "vx_e": 0.25,
        "vx_tmp": 0.25,
        "dth": 0.5,
        "dx": 0.02
    }
    
    fix4 = {
        "vx_e": -0.25,
        "vx_tmp": -0.25,
        "dth": 0.5,
        "dx": 0.02
    }
    
    fix_list = [fix1, fix2, fix3, fix4]
    for fix in fix_list:
        
        i_d, lx = dep_search_1d(
            i, **fix
        )
            
        x_arrival = arrival_point_coordinate(
            i_d, lx,  **fix
        )
                
        assert x_arrival == xcr

if __name__ == "__main__":
    test_dep_search_1d()
