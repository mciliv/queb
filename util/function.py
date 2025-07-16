import inspect

def keys(func: callable):
    params = inspect.signature(func).parameters
    return list(params.keys())

