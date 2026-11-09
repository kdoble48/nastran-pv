"""labels.unique_labels: same-named models told apart by their folders (scene layers rule R6)."""
from nastran_pv.labels import unique_labels


def test_a_unique_file_name_is_the_label():
    assert unique_labels(["D:/runs/a/wing.bdf", "D:/runs/b/fuselage.bdf"]) == ["wing.bdf", "fuselage.bdf"]


def test_the_same_file_name_takes_its_folder():
    paths = [r"D:\ws\tests\data\same_name\a\model.bdf", "D:/ws/tests/data/same_name/b/model.bdf"]
    assert unique_labels(paths) == ["a/model.bdf", "b/model.bdf"]
    ucrm = ["X:/diff_example/A/run_uCRM-9_modal_coarse.bdf", "X:/diff_example/B/run_uCRM-9_modal_coarse.bdf"]
    assert unique_labels(ucrm) == ["A/run_uCRM-9_modal_coarse.bdf", "B/run_uCRM-9_modal_coarse.bdf"]


def test_as_many_folders_as_it_takes_for_every_path_of_that_name():
    paths = ["/r/x/a/m.bdf", "/r/y/a/m.bdf", "/r/z/b/m.bdf", "/r/other.bdf"]
    assert unique_labels(paths) == ["x/a/m.bdf", "y/a/m.bdf", "z/b/m.bdf", "other.bdf"]
    assert unique_labels(["a/m.bdf", "/r/x/a/m.bdf"]) == ["a/m.bdf", "x/a/m.bdf"]  # one is the other's tail


def test_equal_paths_keep_equal_labels():
    """The same deck opened twice: the caller adds the proxy id to the later one."""
    assert unique_labels(["/r/a/m.bdf", "/r/a/m.bdf"]) == ["m.bdf", "m.bdf"]
    assert unique_labels(["/r/a/m.bdf", "/r/A/M.bdf", "/r/b/m.bdf"]) == ["a/m.bdf", "A/M.bdf", "b/m.bdf"]


def test_no_paths():
    assert unique_labels([]) == []
