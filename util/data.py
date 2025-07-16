def merge(priority: dict, secondary: dict) -> dict:
    """Merge two dictionaries, with priority dict taking precedence"""
    result = secondary.copy()
    result.update(priority)
    return result
        
