from functools import partial
from hashlib import sha256

from luigi import ExternalTask
from luigi import LocalTarget
from luigi import Parameter
from luigi import Task
from luigi.task import flatten


def creat_salted_version(task):
    """Create a salted id/version for this task and lineage
    """
    Pass



class Requirement:
    def __init__(self, task_class, **params):
        self.task_class = task_class
        self.params = params

    def __get__(self, task, cls):
        # generate clone
        if task is None:
            return type(task)

        return task.clone(
            self.task_class, **self.params)


class Requires:
    """Composition to replace :meth:`luigi.task.Task.requires`

    Example::

        class MyTask(Task):
            # Replace task.requires()
            requires = Requires()
            other = Requirement(OtherTask)

            def run(self):
                # Convenient access here...
                with self.other.output().open('r') as f:
                    ...

        >>> MyTask().requires()
        {'other': OtherTask()}

    """

    def __get__(self, task, cls):
        if task is None:
            return self
        # Bind self/task in a closure
        return partial(self.__call__, task)

    def __call__(self, task):
        """Returns the requirements of a task

        Assumes the task class has :class:`.Requirement` descriptors, which
        can clone the appropriate dependences from the task instance.

        :returns: requirements compatible with `task.requires()`
        :rtype: dict
        """
        # Search task.__class__ for Requirement instances
        # return
        attribute_dic = {k: getattr(task, k) for k in dir(task)}
        return {k: v for k, v in attribute_dic.items() if isinstance(v, Task)}


class TargetOutput:
    def __init__(self, file_pattern='{task.__class__.__name__}',
                 ext='.txt', target_class=LocalTarget, **target_kwargs):
        self.target_class = target_class
        self.target_kwargs = target_kwargs
        self.file_pattern = file_pattern
        self.ext = ext

    def __get__(self, task, cls):
        return partial(self.__call__, task)

    def __call__(self, task):
        # Determine the path etc here
        file_name = self.file_pattern + self.ext
        return self.target_class(file_name.format(task=task))



