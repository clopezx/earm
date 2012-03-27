def write_dict(kd):
    """ this function assumes that model.parameters already has the right values assigned
    into the model... (which gets automagically propagated to kd). Eventually this will
    be expanded to write to a file.
    """
    # FIXME: this just dumps a list from the OrderedDict
    #        a "prettier" representation would be nice...
    # might need to fix the spacing in the new dict by hand... fortunately this is only done once in a while
    for i in kd:
        print "(",repr(i),","
        print repr(kd[i]), "),"

        
                           
                           
