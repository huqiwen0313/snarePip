from dask import delayed
from dask.bytes.core import get_fs_token_paths
from dask.dataframe import read_csv
from dask.dataframe import read_parquet
from dask.dataframe import to_csv
from dask.dataframe import to_parquet
from luigi import Target
from luigi.task import logger as luigi_logger

from csci_utils.io.hash_str import aws_keys

FLAG = "_SUCCESS"  # Task success flag, c.f. Spark, hadoop, etc


@delayed(pure=True)
def touch(path, storage_options=None, _dep=None):
    fs, token, paths = get_fs_token_paths(path, storage_options=storage_options)
    with fs.open(path, mode="wb"):
        pass


class BaseDaskTarget(Target):
    """Base target for dask collections

    The class provides reading and writing mechanisms to any filesystem
    supported by dask, as well as a standardized method of flagging a task as
    successfully completed (or checking it).

    These targets can be used even if Dask is not desired, as a way of handling
    paths and success flags, since the file structure (parquet/csvs/jsons etc)
    are identical to other Hadoop and alike systems.
    """

    def __init__(self, path, glob=None, flag=FLAG, storage_options=None):
        """

        :param str path: Directory the collection is stored in.  May be remote
            with the appropriate dask backend (ie s3fs)

        :param str flag: file under directory which indicates a successful
            completion, like a hadoop '_SUCCESS' flag.  If no flag is written,
            or another file cannot be used as a proxy, then the task can only
            assume a successful completion based on a file glob being
            representative, ie the directory exists or at least 1 file matching
            the glob. You must set the flag to '' or None to allow this fallback behavior.

        :param str glob: optional glob for files when reading only (or as a
            proxy for completeness) or to set the name for files to write for
            output types which require such.

        :param dict storage_options: used to create dask filesystem or passed on read/write
        """

        self.glob = glob
        self.path = path
        self.storage_options = storage_options or {}
        self.flag = flag

        #if not path.endswith(self._get_sep()):
        #    raise ValueError("Must be a directory!")

    @property
    def fs(self):
        fs, token, paths = get_fs_token_paths(
            self.path, storage_options=self.storage_options
        )
        return fs

    def _get_sep(self):
        """The path seperator for the fs"""
        try:
            return self.fs.sep  # Set for local files etc
        except AttributeError:
            # Typical for s3, hdfs, etc
            return "/"

    def _join(self, *paths):
        sep = self._get_sep()
        pths = [p.rstrip(sep) for p in paths]
        return sep.join(pths)

    def _exists(self, path):
        """Check if some path or glob exists

        This should not be an implementation, yet the underlying FileSystem
        objects do not share an interface!

        :rtype: bool
        """

        try:
            _e = self.fs.exists
        except AttributeError:
            try:
                return len(self.fs.glob(path)) > 0
            except FileNotFoundError:
                return False
        return _e(path)

    def exists(self):
        # NB: w/o a flag, we cannot tell if a task is partially done or complete
        fs = self.fs
        if self.flag:
            ff = self._join(self.path, self.flag)
            if hasattr(fs, 'exists'):
                # Unfortunately, not every dask fs implemenents this!
                return fs.exists(ff)

        else:
            ff = self._join(self.path, self.glob or "")

        try:
            for _ in fs.glob(ff):
                # If a single file is found, assume exists
                # for loop in case glob is iterator, no need to consume all
                return True

            # Empty list?
            return False

        except FileNotFoundError:
            return False

    def mark_complete(self, _dep=None, compute=True):
        if not self.flag:
            raise RuntimeError("No flag for task, cannot mark complete")

        flagfile = self._join(self.path, self.flag)

        out = touch(flagfile, _dep=_dep)
        if compute:
            out.compute()
            return
        return out

    def augment_options(self, storage_options):
        """Get composite storage options

        :param dict storage_options: for a given call, these will take precedence

        :returns: options with defaults baked in
        :rtype: dict
        """
        base = self.storage_options.copy()
        if storage_options is not None:
            base.update(storage_options)
        return base

    def read_dask(self, storage_options=None, check_complete=True, **kwargs):
        #if check_complete and self.flag and not self.exists():
        #    raise FileNotFoundError("Task not yet run or incomplete")

        return self._read(
            self.get_path_for_read(),
            storage_options=self.augment_options(storage_options),
            **kwargs
        )

    def get_path_for_read(self):
        if self.glob:
            return self._join(self.path, self.glob)
        return self.path.rstrip(self._get_sep())

    def get_path_for_write(self):
        return self.path.rstrip(self._get_sep())

    def write_dask(
        self,
        collection,
        compute=True,
        storage_options=None,
        logger=luigi_logger,
        **kwargs
    ):
        if logger:
            logger.info("Writing dask collection to {}".format(self.path))
        storage_options = self.augment_options(storage_options)
        out = self._write(
            collection,
            self.get_path_for_write(),
            storage_options=storage_options,
            compute=False,
            **kwargs
        )

        if self.flag:
            out = self.mark_complete(_dep=out, compute=False)

        if compute:
            out.compute()
            if logger:
                logger.info(
                    "Successfully wrote to {}, flagging complete".format(self.path)
                )
            return None
        return out

    # Abstract interface: implement for various storage formats

    @classmethod
    def _read(cls, path, **kwargs):
        raise NotImplementedError()

    @classmethod
    def _write(cls, collection, path, **kwargs):
        raise NotImplementedError()


class ParquetTarget(BaseDaskTarget):
    @classmethod
    def _read(cls, path, **kwargs):
        cls.path = path
        cls.kwargs = kwargs
        data = read_parquet(cls.path, **kwargs)
        return data

    @classmethod
    def _write(cls, collection, path, **kwargs):
        cls.collection = collection
        cls.path = path
        cls.collection.to_parquet(cls.path)


class CSVTarget(BaseDaskTarget):
    # read file and convert into dask dataframe
    @classmethod
    def _read(cls, path, **kwargs):
        cls.path = path
        cls.kwargs = kwargs
        data = read_csv(cls.path, **kwargs)
        return data

    @classmethod
    def _write(cls, collection, path, **kwargs):
        cls.collection = collection
        cls.path = path
        cls.collection.to_csv(cls.path)
