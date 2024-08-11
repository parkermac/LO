"""
This is a module of experiment-specific customizations that are applied
within tracker2.

It is not clear to me the best way to implement such specific behaviors in the
code, so this module is an attempt.
"""

def update_TR(nd, TR):
    """
    This was motivated by the Nina Bednarsek experiments where she wanted larvae
    (rockfish?) to stay at one depth for the first month and then go to a different
    depth for the second month.
    """
    # Note that because TR is a dictionary, changing it in the function
    # changes it everywhere. No need to return it.
    
    if ('nina_' in TR['exp_name']) and (TR['sub_tag'] == 'twopart') and (nd >= 30):
        TR['stay'] = 5000 # 0 causes stay depth to be ignored
    

