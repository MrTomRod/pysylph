from . import test_database, test_doctest, test_sketcher, test_profiler, test_twostage


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_database))
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_sketcher))
    suite.addTests(loader.loadTestsFromModule(test_profiler))
    suite.addTests(loader.loadTestsFromModule(test_twostage))
    return suite
