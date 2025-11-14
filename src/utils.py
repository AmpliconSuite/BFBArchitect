import logging

def create_logger(name, log_file):
    """Create a logger"""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    
    logger.handlers.clear()
    
    # Create file handler
    handler = logging.FileHandler(log_file, mode='w')
    handler.setLevel(logging.DEBUG)
    
    # Create formatter
    formatter = logging.Formatter('[%(name)s:%(levelname)s]\t%(message)s')
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    return logger