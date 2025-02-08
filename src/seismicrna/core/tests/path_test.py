import unittest as ut
from shutil import rmtree

from seismicrna.core.path import *


class TestGetSeismicRNASourceDir(ut.TestCase):

    def test_get_seismicrna_source_dir(self):
        seismicrna_source_dir = get_seismicrna_source_dir()
        self.assertEqual(seismicrna_source_dir,
                         sanitize(__file__).parent.parent.parent)


class TestValidateBranches(ut.TestCase):

    def test_validate_branches(self):
        self.assertIsNone(validate_branches({"step1": ""}))
        self.assertIsNone(validate_branches({"step1": "a1"}))
        self.assertIsNone(validate_branches({"step_1": "a1",
                                             "step_2": "b1"}))
        self.assertRaisesRegex(PathTypeError,
                               "branches must be dict, but got list",
                               validate_branches,
                               ["a1"])
        self.assertRaisesRegex(PathTypeError,
                               "txt must be str, but got int",
                               validate_branches,
                               {1: "a1"})
        self.assertRaisesRegex(PathTypeError,
                               "branch must be str, but got int",
                               validate_branches,
                               {"step1": 1})
        self.assertRaisesRegex(PathValueError,
                               "txt cannot be empty string",
                               validate_branches,
                               {"": "a1"})
        self.assertRaisesRegex(PathValueError,
                               r"txt 'step 1' has illegal characters: \[' '\]",
                               validate_branches,
                               {"step 1": "a1"})
        self.assertRaisesRegex(PathValueError,
                               r"branch 'a_1' has illegal characters: \['_'\]",
                               validate_branches,
                               {"step1": "a_1"})

    def test_validate_branches_flat(self):
        self.assertIsNone(validate_branches_flat([]))
        self.assertIsNone(validate_branches_flat(["a1"]))
        self.assertIsNone(validate_branches_flat(["a1", "b1"]))
        self.assertRaisesRegex(PathTypeError,
                               "branches_flat must be list, but got dict",
                               validate_branches_flat,
                               {"step1": "a1"})
        self.assertRaisesRegex(PathTypeError,
                               "branch must be str, but got int",
                               validate_branches_flat,
                               [1])
        self.assertRaisesRegex(PathValueError,
                               "branch cannot be the empty string",
                               validate_branches_flat,
                               [""])
        self.assertRaisesRegex(PathValueError,
                               r"branch 'a 1' has illegal characters: \[' '\]",
                               validate_branches_flat,
                               ["a 1"])
        self.assertRaisesRegex(PathValueError,
                               r"branch 'a_1' has illegal characters: \['_'\]",
                               validate_branches_flat,
                               ["a_1"])


class TestAddBranch(ut.TestCase):

    def test_add_branches(self):
        self.assertDictEqual(add_branch("step1", "", {}),
                             {"step1": ""})
        self.assertDictEqual(add_branch("step2", "", {"step1": "a1"}),
                             {"step1": "a1", "step2": ""})
        self.assertDictEqual(add_branch("step1", "b1", {}),
                             {"step1": "b1"})
        self.assertDictEqual(add_branch("step2", "b1", {"step1": "a1"}),
                             {"step1": "a1",
                              "step2": "b1"})
        self.assertDictEqual(add_branch("step3", "b1", {"step1": "a1",
                                                        "step2": "a2"}),
                             {"step1": "a1",
                              "step2": "a2",
                              "step3": "b1"})
        self.assertRaisesRegex(PathValueError,
                               "txt cannot be empty string",
                               add_branch,
                               "", "a1", {})
        self.assertRaisesRegex(PathValueError,
                               r"branch 'a_1' has illegal characters: \['_'\]",
                               add_branch,
                               "step1", "a_1", {})
        self.assertRaisesRegex(PathValueError,
                               "Duplicate step: 'step1'",
                               add_branch,
                               "step1", "a2", {"step1": "a1"})


class TestBranchesField(ut.TestCase):

    def test_branches_validate(self):
        self.assertIsNone(BranchesField.validate([]))
        self.assertIsNone(BranchesField.validate(["branch1"]))
        self.assertIsNone(BranchesField.validate(["branch1", "branch2"]))
        self.assertRaisesRegex(PathTypeError,
                               "branches_flat must be list, but got str",
                               BranchesField.validate,
                               "branch1")
        self.assertRaisesRegex(PathTypeError,
                               "branch must be str, but got int",
                               BranchesField.validate,
                               [1])
        self.assertRaisesRegex(
            PathValueError,
            "branch cannot be the empty string",
            BranchesField.validate,
            [""]
        )
        self.assertRaisesRegex(
            PathValueError,
            r"branch 'branch 1' has illegal characters: \[' '\]",
            BranchesField.validate,
            ["branch 1"]
        )
        self.assertRaisesRegex(
            PathValueError,
            r"branch 'branch_1' has illegal characters: \['_'\]",
            BranchesField.validate,
            ["branch_1"]
        )

    def test_branches_build(self):
        self.assertEqual(BranchesField.build([]),
                         "")
        self.assertEqual(BranchesField.build(["branch1"]),
                         "_branch1")
        self.assertEqual(BranchesField.build(["branch1"]),
                         "_branch1")
        self.assertEqual(BranchesField.build(["branch1", "branch2"]),
                         "_branch1_branch2")
        self.assertRaisesRegex(PathTypeError,
                               "branches_flat must be list, but got str",
                               BranchesField.build,
                               "branch1")
        self.assertRaisesRegex(PathValueError,
                               "branch cannot be the empty string",
                               BranchesField.build,
                               [""])

    def test_branches_parse(self):
        self.assertListEqual(BranchesField.parse(""),
                             [])
        self.assertListEqual(BranchesField.parse("branch-1"),
                             ["branch-1"])
        self.assertListEqual(BranchesField.parse("_branch-1"),
                             ["branch-1"])
        self.assertListEqual(BranchesField.parse("branch-1_"),
                             ["branch-1"])
        self.assertListEqual(BranchesField.parse("_branch-1_"),
                             ["branch-1"])
        self.assertListEqual(BranchesField.parse("_branch1_branch2"),
                             ["branch1", "branch2"])
        self.assertListEqual(BranchesField.parse("_branch1__branch2"),
                             ["branch1", "branch2"])
        self.assertListEqual(BranchesField.parse("_branch1__branch2_"),
                             ["branch1", "branch2"])
        self.assertListEqual(BranchesField.parse("_branch1___branch2_"),
                             ["branch1", "branch2"])
        self.assertRaisesRegex(PathTypeError,
                               "text to parse must be str, but got list",
                               BranchesField.parse,
                               ["branch1"])
        self.assertRaisesRegex(
            PathValueError,
            r"branch 'branch 2' has illegal characters: \[' '\]",
            BranchesField.parse,
            "_branch1_branch 2"
        )


class TestCmdSeg(ut.TestCase):

    def test_cmdseg_build_branches_list(self):
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: []}),
                         "align")
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: ["branch1"]}),
                         "align_branch1")
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: ["branch1", "branch2"]}),
                         "align_branch1_branch2")
        self.assertRaisesRegex(PathValueError,
                               "branch cannot be the empty string",
                               CmdSeg.build,
                               {CMD: "align",
                                BRANCHES: ["branch1", ""]})
        self.assertRaisesRegex(PathValueError,
                               "Invalid option 'malign'",
                               CmdSeg.build,
                               {CMD: "malign",
                                BRANCHES: ["branch1"]})

    def test_cmdseg_build_branches_dict(self):
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: {}}),
                         "align")
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: {"step1": ""}}),
                         "align")
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: {"step1": "branch1"}}),
                         "align_branch1")
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: {"step1": "",
                                                  "step2": "branch1"}}),
                         "align_branch1")
        self.assertEqual(CmdSeg.build({CMD: "align",
                                       BRANCHES: {"step1": "branch1",
                                                  "step2": "branch2"}}),
                         "align_branch1_branch2")
        self.assertRaisesRegex(PathValueError,
                               "Invalid option 'malign'",
                               CmdSeg.build,
                               {CMD: "malign",
                                BRANCHES: {"step1": "branch1"}})

    def test_cmdseg_parse(self):
        self.assertEqual(CmdSeg.parse("align"),
                         {CMD: "align", BRANCHES: []})
        self.assertEqual(CmdSeg.parse("align_branch1"),
                         {CMD: "align", BRANCHES: ["branch1"]})
        self.assertEqual(CmdSeg.parse("align__branch1"),
                         {CMD: "align", BRANCHES: ["branch1"]})
        self.assertEqual(CmdSeg.parse("align_branch1_"),
                         {CMD: "align", BRANCHES: ["branch1"]})
        self.assertEqual(CmdSeg.parse("align_branch1_branch2"),
                         {CMD: "align", BRANCHES: ["branch1", "branch2"]})
        self.assertEqual(CmdSeg.parse("align__branch1_branch2"),
                         {CMD: "align", BRANCHES: ["branch1", "branch2"]})
        self.assertEqual(CmdSeg.parse("align_branch1__branch2"),
                         {CMD: "align", BRANCHES: ["branch1", "branch2"]})
        self.assertRaisesRegex(PathValueError,
                               "Could not parse fields in text ''",
                               CmdSeg.parse,
                               "")
        self.assertRaisesRegex(PathValueError,
                               "Could not parse fields in text '_align'",
                               CmdSeg.parse,
                               "_align")
        self.assertRaisesRegex(PathValueError,
                               "Invalid option 'alight'",
                               CmdSeg.parse,
                               "alight")


class TestSymlinkIfNeeded(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._tmp_dir = None

    def setUp(self):
        self._tmp_dir = randdir()

    def tearDown(self):
        rmtree(self._tmp_dir)
        self._tmp_dir = None

    def test_target_not_exist(self):
        link = self._tmp_dir.joinpath("link")
        target = self._tmp_dir.joinpath("target")
        for link_exists in range(2):
            self.assertRaisesRegex(FileNotFoundError,
                                   str(target),
                                   symlink_if_needed,
                                   link,
                                   target)
            link.mkdir(exist_ok=link_exists)
            self.assertTrue(link.is_dir())

    def test_link_valid(self):
        link = self._tmp_dir.joinpath("link")
        target = self._tmp_dir.joinpath("target")
        target.mkdir()
        for _ in range(2):
            symlink_if_needed(link, target)
            self.assertTrue(link.is_symlink())
            self.assertTrue(link.readlink() == target)

    def test_link_not_symlink(self):
        link = self._tmp_dir.joinpath("link")
        link.mkdir()
        target = self._tmp_dir.joinpath("target")
        target.mkdir()
        self.assertRaisesRegex(OSError,
                               f"{link} is not a symbolic link",
                               symlink_if_needed,
                               link,
                               target)

    def test_link_wrong_symlink(self):
        link = self._tmp_dir.joinpath("link")
        target = self._tmp_dir.joinpath("target")
        target.mkdir()
        target2 = self._tmp_dir.joinpath("target2")
        link.symlink_to(target2)
        self.assertRaisesRegex(
            OSError,
            f"{link} is a symbolic link to {target2}, not to {target}",
            symlink_if_needed,
            link,
            target
        )


if __name__ == "__main__":
    ut.main(verbosity=2)
