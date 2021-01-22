''' Custom masks

'''


def topo1(hy_obj):
    ir = hy_obj.get_wave(850)
    red = hy_obj.get_wave(660)
    ndvi = (ir-red)/(ir+red)
    cosine_i = hy_obj.cosine_i()
    slope = hy_obj.get_anc('slope')
    mask = (ndvi > 0.4) & (ir != hy_obj.no_data) & (cosine_i > .12) & (slope > .03)
    return mask

def brdf1(hy_obj):
    ir = hy_obj.get_wave(850)
    red = hy_obj.get_wave(660)
    ndvi = (ir-red)/(ir+red)
    mask = (ir != hy_obj.no_data)
    return mask



mask_dict = {'topo1' : topo1,
              'brdf1' : brdf1}
              

