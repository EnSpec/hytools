"""Miscellaneous functions

"""
from itertools import tee

def progbar(curr, total, full_progbar = 100):
    '''Display progress bar.

    Gist from:

    https://gist.github.com/marzukr/3ca9e0a1b5881597ce0bcb7fb0adc549

    Args:
        curr (int, float): Current task level.
        total (int, float): Task level at completion.
        full_progbar (TYPE): Defaults to 100.

    Returns:
        None.

    '''
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')


def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def set_brdf(hy_obj,brdf_dict):
    hy_obj.brdf = brdf_dict

def update_brdf(hy_obj,args):
    hy_obj.brdf[args['key']] = args['value']