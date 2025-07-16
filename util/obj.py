import util.dictionary as d


def eq_attributes(first: type, second: type, ignore=[]):
    return dictionary(first, ignore) == dictionary(second, ignore)

def dictionary(obj: type, ignore=[]):
    return d.remove(vars(obj), ignore)
