from astropy import logger
import logging

__all__ = ['pirpyfromat','log','create_log']

pirpyfromat = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    
class PirpyLogger(logger.AstropyLogger):
    def enable_log_to_file(self, filename, filter_level=None, filter_origin=None,
                           format=None):
        self._fh = logging.FileHandler(filename)
        if filter_level is not None:
            self._fh.setLevel(filter_level)
        if filter_origin is not None:
            self._fh.addFilter(FilterOrigin(filter_origin))
        if format is not None:
            f = logging.Formatter(format)
            self._fh.setFormatter(f)
        self.addHandler(self._fh)
        
    def disable_log_to_file(self):
        self._fh.close()
        self.removeHandler(self._fh)

def create_log(name):
    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(PirpyLogger)
    try:
        log = logging.getLogger(name)
        log._set_defaults()
    finally:
        logging.setLoggerClass(orig_logger_cls)

    return log

log = create_log('pirpy')
log.setLevel('INFO')
