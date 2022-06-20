# taken from https://xsnippet.org/359377/

# TODO: modification https://stackoverflow.com/questions/20333674/pycharm-logging-output-colours
#  but I have already handled the case (I guess..)

import sys
import logging

CEND = '\33[0m'
CBOLD = '\33[1m'
CITALIC = '\33[3m'
CURL = '\33[4m'
CBLINK = '\33[5m'
CBLINK2 = '\33[6m'
CSELECTED = '\33[7m'

CBLACK = '\33[30m'
CRED = '\33[31m'
CGREEN = '\33[32m'
CYELLOW = '\33[33m'
CBLUE = '\33[34m'
CVIOLET = '\33[35m'
CBEIGE = '\33[36m'
CWHITE = '\33[37m'

CBLACKBG = '\33[40m'
CREDBG = '\33[41m'
CGREENBG = '\33[42m'
CYELLOWBG = '\33[43m'
CBLUEBG = '\33[44m'
CVIOLETBG = '\33[45m'
CBEIGEBG = '\33[46m'
CWHITEBG = '\33[47m'

CGREY = '\33[90m'
CRED2 = '\33[91m'
CGREEN2 = '\33[92m'
CYELLOW2 = '\33[93m'
CBLUE2 = '\33[94m'
CVIOLET2 = '\33[95m'
CBEIGE2 = '\33[96m'
CWHITE2 = '\33[97m'

CGREYBG = '\33[100m'
CREDBG2 = '\33[101m'
CGREENBG2 = '\33[102m'
CYELLOWBG2 = '\33[103m'
CBLUEBG2 = '\33[104m'
CVIOLETBG2 = '\33[105m'
CBEIGEBG2 = '\33[106m'
CWHITEBG2 = '\33[107m'


class _AnsiColorizer(object):
    """
    A colorizer is an object that loosely wraps around a stream, allowing
    callers to write text to the stream in a particular color.

    Colorizer classes must implement C{supported()} and C{write(text, color)}.
    """
    _colors = dict(black=30, red=31, green=32, yellow=33,
                   blue=34, magenta=35, cyan=36, white=37)

    # _colors = dict(black=30, red=31, green=32, yellow=33,
    #                blue=34, magenta=35, cyan=36, white=37)

    def __init__(self, stream):
        self.stream = stream

    @classmethod
    def supported(cls, stream=sys.stdout):
        """
        A class method that returns True if the current platform supports
        coloring terminal output using this method. Returns False otherwise.
        """
        if not stream.isatty():
            return False  # auto color only on TTYs
        try:
            import curses
        except ImportError:
            return False
        else:
            try:
                try:
                    return curses.tigetnum("colors") > 2
                except curses.error:
                    curses.setupterm()
                    return curses.tigetnum("colors") > 2
            except:
                raise
                # guess false in case of error
                return False

    def write(self, text, color):
        """
        Write the given text to the stream in the given color.

        @param text: Text to be written to the stream.

        @param color: A string label for a color. e.g. 'red', 'white'.
        """
        # color = self._colors[color]
        # self.stream.write('\x1b[%s;1m%s\x1b[0m' % (color, text))
        # style = 0
        # formatting = ';'.join([str(style), str(color), str(40)])
        # self.stream.write('\x1b[%s;1m%s\x1b[0m' % (color, text))
        self.stream.write(f'{color}{text}\x1b[0m')
        # for style in range(8):
        #     formatting = ';'.join([str(style), str(color), str(40)])
        #     print('\x1b[%sm %s \x1b[0m' % (formatting, text))


class ColorHandler(logging.StreamHandler):
    def __init__(self, stream=sys.stderr):
        super(ColorHandler, self).__init__(_AnsiColorizer(stream))

    def emit(self, record):
        # msg_colors = {
        #     logging.DEBUG: "green",
        #     logging.INFO: "blue",
        #     logging.WARNING: "yellow",
        #     logging.ERROR: "red"
        # }
        msg_colors = {
            logging.DEBUG: "\33[90m",
            logging.INFO: "\33[94m",
            logging.WARNING: "\33[33m",
            logging.ERROR: "\33[31m",
            logging.CRITICAL: "\33[31m"
        }

        # color = msg_colors.get(record.levelno, "blue")
        self.stream.write(record.msg + "\n", msg_colors[record.levelno])

        # logging.getLogger().setLevel(logging.DEBUG)
        # logging.getLogger().addHandler(ColorHandler())
        #
        # if __name__ == "__main__":
        #     logging.debug("Some debugging output")
        #     logging.info("Some info output")
        #     logging.error("Some error output")
        #     logging.warning("Some warning output")
