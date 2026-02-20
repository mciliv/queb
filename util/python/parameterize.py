def parameterize(arg_list):
    def decorator(test_func):
        def wrapper(self):
            for args in arg_list:
                if isinstance(args, list):
                    test_func(self, args)
                else:
                    test_func(self, [args])
        return wrapper
    return decorator
