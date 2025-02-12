import unittest as ut

from seismicrna.core.validate import (require_issubclass,
                                      require_isinstance,
                                      require_isin,
                                      require_equal,
                                      require_atleast,
                                      require_atmost)


class TestRequireIsSubclass(ut.TestCase):

    def test_is_subclass(self):
        classes = OSError
        self.assertIsNone(require_issubclass("xyz",
                                             FileExistsError,
                                             classes))

    def test_is_subclasses(self):
        classes = OSError, int
        self.assertIsNone(require_issubclass("xyz",
                                             FileExistsError,
                                             classes))
        self.assertIsNone(require_issubclass("xyz",
                                             int,
                                             classes))

    def test_not_issubclass(self):
        classes = OSError
        self.assertRaisesRegex(ValueError,
                               (f"xyz must be a subclass of {classes}, "
                                f"but got {ZeroDivisionError}"),
                               require_issubclass,
                               "xyz", ZeroDivisionError, classes)

    def test_not_issubclasses(self):
        classes = OSError, int
        self.assertRaisesRegex(
            ValueError,
            (f"xyz must be a subclass of {classes}, "
             f"but got {float}").replace("(", "[(]").replace(")", "[)]"),
            require_issubclass,
            "xyz", float, classes
        )

    def test_custom_valueerror(self):
        class MyCustomError(ValueError):
            pass

        classes = OSError
        self.assertRaisesRegex(
            MyCustomError,
            (f"xyz must be a subclass of {classes}, "
             f"but got {ZeroDivisionError}"),
            require_issubclass,
            "xyz", ZeroDivisionError, classes, MyCustomError
        )

    def test_custom_not_valueerror(self):
        class MyCustomError(TypeError):
            pass

        classes = OSError
        self.assertRaisesRegex(
            ValueError,
            (f"error_type must be a subclass of {ValueError}, "
             f"but got {MyCustomError}"),
            require_issubclass,
            "xyz", ZeroDivisionError, classes, MyCustomError
        )


class TestRequireIsInstance(ut.TestCase):

    def test_is_instance(self):
        classes = int
        self.assertIsNone(require_isinstance("xyz",
                                             12,
                                             classes))

    def test_is_instances(self):
        classes = int, float
        self.assertIsNone(require_isinstance("xyz",
                                             12,
                                             classes))
        self.assertIsNone(require_isinstance("xyz",
                                             12.3,
                                             classes))

    def test_not_isinstance(self):
        classes = int
        self.assertRaisesRegex(TypeError,
                               (f"xyz must be an instance of {classes}, "
                                f"but got 12.3 of type {float}"),
                               require_isinstance,
                               "xyz", 12.3, classes)

    def test_not_isinstances(self):
        classes = int, float
        self.assertRaisesRegex(
            TypeError,
            (f"xyz must be an instance of {classes}, but got '12.3' "
             f"of type {str}").replace("(", "[(]").replace(")", "[)]"),
            require_isinstance,
            "xyz", "12.3", classes
        )

    def test_custom_typeerror(self):
        class MyCustomError(TypeError):
            pass

        classes = int
        self.assertRaisesRegex(
            MyCustomError,
            (f"xyz must be an instance of {classes}, "
             f"but got 12.3 of type {float}"),
            require_isinstance,
            "xyz", 12.3, classes, MyCustomError
        )

    def test_custom_not_valueerror(self):
        class MyCustomError(ValueError):
            pass

        classes = int
        self.assertRaisesRegex(
            ValueError,
            (f"error_type must be a subclass of {TypeError}, "
             f"but got {MyCustomError}"),
            require_isinstance,
            "xyz", 12.3, classes, MyCustomError
        )


class TestRequireIsIn(ut.TestCase):

    def test_isin(self):
        self.assertIsNone(require_isin("xyz", "b", "abc"))

    def test_not_isin_named(self):
        self.assertRaisesRegex(ValueError,
                               ("xyz must be in letters, "
                                "but got xyz='d' and letters='abc'"),
                               require_isin,
                               "xyz", "d", "abc", "letters")

    def test_not_isin_unnamed(self):
        self.assertRaisesRegex(ValueError,
                               "xyz must be in 'abc', but got 'd'",
                               require_isin,
                               "xyz", "d", "abc")

    def test_custom_exception(self):
        class MyCustomError(Exception):
            pass

        self.assertRaisesRegex(MyCustomError,
                               "xyz must be in 'abc', but got 'd'",
                               require_isin,
                               "xyz", "d", "abc", error_type=MyCustomError)

    def test_custom_not_exception(self):
        class MyCustomError(BaseException):
            pass

        self.assertRaisesRegex(ValueError,
                               ("error_type must be a subclass of "
                                f"{Exception}, but got {repr(MyCustomError)}"),
                               require_isin,
                               "xyz", "d", "abc", error_type=MyCustomError)


class TestRequireEqual(ut.TestCase):

    def test_equal_int(self):
        self.assertIsNone(require_equal("a", 0, 0))

    def test_equal_str(self):
        self.assertIsNone(require_equal("a", "A", "A"))

    def test_equal_wrong_type_value(self):
        self.assertRaisesRegex(TypeError,
                               (f"a must be an instance of {int}, "
                                f"but got '0' of type {str}"),
                               require_equal,
                               "a", "0", 1, classes=int)

    def test_equal_wrong_type_other_unnamed(self):
        self.assertRaisesRegex(TypeError,
                               (f"other must be an instance of {int}, "
                                f"but got '1' of type {str}"),
                               require_equal,
                               "a", 0, "1", classes=int)

    def test_equal_wrong_type_other_named(self):
        self.assertRaisesRegex(TypeError,
                               (f"b must be an instance of {int}, "
                                f"but got '1' of type {str}"),
                               require_equal,
                               "a", 0, "1", other_name="b", classes=int)

    def test_not_equal_int_unnamed(self):
        self.assertRaisesRegex(ValueError,
                               "Must have a = 1, but got 0",
                               require_equal,
                               "a", 0, 1)

    def test_not_equal_str_unnamed(self):
        self.assertRaisesRegex(ValueError,
                               "Must have a = 'A', but got 'B'",
                               require_equal,
                               "a", "B", "A")

    def test_not_equal_int_named(self):
        self.assertRaisesRegex(ValueError,
                               "Must have a = b, but got a=0 and b=1",
                               require_equal,
                               "a", 0, 1, other_name="b")

    def test_not_equal_str_named(self):
        self.assertRaisesRegex(ValueError,
                               "Must have a = b, but got a='A' and b='B'",
                               require_equal,
                               "a", "A", "B", other_name="b")

    def test_custom_valueerror(self):
        class MyCustomError(ValueError):
            pass

        self.assertRaisesRegex(MyCustomError,
                               "Must have a = 1, but got 0",
                               require_equal,
                               "a", 0, 1, error_type=MyCustomError)

    def test_custom_not_valueerror(self):
        class MyCustomError(TypeError):
            pass

        self.assertRaisesRegex(ValueError,
                               ("error_type must be a subclass of "
                                f"{ValueError}, but got {repr(MyCustomError)}"),
                               require_equal,
                               "a", 0, 1, error_type=MyCustomError)


class TestRequireAtLeast(ut.TestCase):

    def test_atleast(self):
        self.assertIsNone(require_atleast("a", 0, 0))
        self.assertIsNone(require_atleast("a", 1, 0))

    def test_not_atmost(self):
        self.assertRaisesRegex(ValueError,
                               "Must have a ≥ 1, but got 0",
                               require_atleast,
                               "a", 0, 1)

    def test_atleast_wrong_type_minimum_unnamed(self):
        self.assertRaisesRegex(TypeError,
                               (f"minimum must be an instance of {int}, "
                                f"but got '0' of type {str}"),
                               require_atleast,
                               "a", 1, "0", classes=int)


class TestRequireAtMost(ut.TestCase):

    def test_atmost(self):
        self.assertIsNone(require_atmost("a", 0, 0))
        self.assertIsNone(require_atmost("a", 0, 1))

    def test_not_atmost(self):
        self.assertRaisesRegex(ValueError,
                               "Must have a ≤ 0, but got 1",
                               require_atmost,
                               "a", 1, 0)

    def test_atmost_wrong_type_maximum_unnamed(self):
        self.assertRaisesRegex(TypeError,
                               (f"maximum must be an instance of {int}, "
                                f"but got '1' of type {str}"),
                               require_atmost,
                               "a", 0, "1", classes=int)
