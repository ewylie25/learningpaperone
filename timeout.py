## {{{ http://code.activestate.com/recipes/577853/ (r1)
"""
Code to timeout with processes.

>>> @timeout(.5)
... def sleep(x):
...     print "ABOUT TO SLEEP {0} SECONDS".format(x)
...     time.sleep(x)
...     return x

>>> sleep(1)
Traceback (most recent call last):
   ...
TimeoutException: timed out after 0 seconds

>>> sleep(.2)
0.2

>>> @timeout(.5)
... def exc():
...     raise Exception('Houston we have problems!')

>>> exc()
Traceback (most recent call last):
   ...
Exception: Houston we have problems!

"""
import multiprocessing
import time
import logging
import traceback
logger = multiprocessing.log_to_stderr()
logger.setLevel(logging.CRITICAL)


class TimeoutException(Exception):
    pass


class RunableProcessing(multiprocessing.Process):
    def __init__(self, func, *args, **kwargs):
        self.queue = multiprocessing.Queue(maxsize=1)
        args = (func,) + args
        multiprocessing.Process.__init__(self, target=self.run_func, args=args, kwargs=kwargs)

    def run_func(self, func, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
            self.queue.put((True, result))
        except Exception as e:
            if "Boost.Python" in str(e.__class__):
                e = Exception(str(e.__class__) + " " +
                              e.message + traceback.format_exc(e))
            self.queue.put((False, e))

    def done(self):
        return self.queue.full()

    def result(self):
        return self.queue.get()


def timeout(seconds, force_kill=True):
    def wrapper(function):
        def inner(*args, **kwargs):
            now = time.time()
            proc = RunableProcessing(function, *args, **kwargs)
            proc.start()
            proc.join(seconds)
            if proc.is_alive():
                if force_kill:
                    proc.terminate()
                runtime = int(time.time() - now)
                raise TimeoutException('Function {0} timed out after {1} seconds\n args: {2}, kwargs: {3} '.format(function.func_name, runtime, args, kwargs))
            assert proc.done()
            success, result = proc.result()
            if success:
                return result
            else:
                raise result
        return inner
    return wrapper


if __name__ == '__main__':
    import doctest
    doctest.testmod()
## end of http://code.activestate.com/recipes/577853/ }}}
