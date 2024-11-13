import signal


class timeout:
    # This is NOT thread safe (only one clocking timeout per thread)
    # https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish
    # slightly modified

    def __init__(self, seconds: int | None = 1, error: Exception = TimeoutError('')):
        self.seconds = seconds
        self.error = error

    def handle_timeout(self, signum, frame):
        raise self.error

    def __enter__(self):
        if self.seconds is not None and self.seconds > 0:
            signal.signal(signal.SIGALRM, self.handle_timeout)
            signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)