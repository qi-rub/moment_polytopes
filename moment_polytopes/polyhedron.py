from __future__ import absolute_import, print_function
import os, subprocess, tempfile, errno, logging
from sage.all import vector, QQ, Polyhedron

__all__ = ["HRepr", "VRepr", "LrsError"]

logger = logging.getLogger(__name__)


class LrsError(Exception):
    """An error that occurred while invoking ``lrslib``."""

    def __init__(self, cmd, message):
        super(Exception, self).__init__('Error while invoking "%s": %s' %
                                        (cmd, message))
        #: Command that was invoked.
        self.cmd = cmd

        #: Error message.
        self.message = message


def _call_lrs(cmd, stdin):
    """Call an command line utility from the lrs suite. Use ``sage -i lrslib`` to install."""

    # NOTE: :func:`subprocess.communicate` runs into trouble for large outputs, hence we do the following

    # generate temporary paths
    def temp_path(name):
        return os.path.join(tempfile.gettempdir(),
                            'justamoment-lrs-%s-%s.txt' % (os.getpid(), name))

    path_stdin = temp_path('stdin')
    path_stdout = temp_path('stdout')
    path_stderr = temp_path('stderr')
    open(path_stdin, 'w').write(stdin)

    # run lrs command
    logger.debug('Running "%s" (%s)', cmd, path_stdin)
    try:
        p = subprocess.Popen(
            cmd,
            stdin=open(path_stdin, 'r'),
            stdout=open(path_stdout, 'w'),
            stderr=open(path_stderr, 'w'))
    except OSError as e:
        if e.errno == errno.ENOENT:
            raise LrsError(cmd, "Not found.")
        else:
            raise
    p.wait()
    logger.debug('"%s" done (%s)', cmd, path_stdout)
    stdout = open(path_stdout, 'r').read().strip()
    stderr = open(path_stderr, 'r').read().strip()
    if stderr:
        raise LrsError(cmd, stderr)

    # cleanup if no error occurred
    os.unlink(path_stdin)
    os.unlink(path_stdout)
    os.unlink(path_stderr)
    return stdout


class HRepr(object):
    """Light-weight container that stores a rational polyhedron in H-representation.

    Inequalities

    .. math::

        H_1 x_1 + \dots + H_d x_d \geq c

    are represented by pairs :math:`(H, c)`, where :math:`H` is a :class:`vector` and :math:`c` an integer (and likewise for equations).

    A H-representation without any inequalities and equations is assumed to be the full ambient space.
    Inequalities and equations are sorted lexicographically.

    :param ieqs: inequalities.
    :param eqns: equations.
    :param ambient_dim: ambient dimension. Needs to be specified only if both ``ieqs`` and ``eqns`` are empty.
    :raises ValueError: neither ``ieqs``, ``eqns``, nor ``ambient_dim`` were specified.

    .. automethod:: __contains__
    .. automethod:: __eq__
    .. automethod:: __and__
    """

    def __init__(self, ieqs=[], eqns=[], ambient_dim=None):
        # convert H to sage vectors/scalars, and sort lexicographically
        def convert((H, c)):
            return (vector(QQ, H), QQ(c))

        #: The inequalities.
        self.ieqs = sorted(convert(ieq) for ieq in ieqs)

        #: The equations.
        self.eqns = sorted(convert(eqn) for eqn in eqns)

        # determine ambient dimension
        if not (ambient_dim or (ieqs or eqns)):
            raise ValueError(
                'Need to specify ambient dimension if ieqs and eqns are empty.')
        if not ambient_dim:
            H, _ = next(iter(self.ieqs if self.ieqs else self.eqns))
            ambient_dim = len(H)

        #: The ambient dimension.
        self.ambient_dim = ambient_dim

        # check dimension
        assert all(len(H) == ambient_dim for (H, _) in self.ieqs)
        assert all(len(H) == ambient_dim for (H, _) in self.eqns)

    def __str__(self):
        return '<HRepr object at 0x%016X in dimension %d defined by %d inequalities and %d equations>' % (
            id(self), self.ambient_dim, len(self.ieqs), len(self.eqns))

    def __eq__(self, rhs):
        """Compare H-representations.

        *Warning:* Equal polyhedra can have different H-representation!

        :rtype: bool"""
        return self.ambient_dim == rhs.ambient_dim and self.ieqs == rhs.ieqs and self.eqns == rhs.eqns

    def __contains__(self, pt):
        """Check if given point ``pt`` is contained in the polyhedron.

        :rtype: bool
        """
        pt = vector(pt)
        return all(pt.dot_product(H) >= c for (H, c) in self.ieqs) \
            and all(pt.dot_product(H) == c for (H, c) in self.eqns)

    def vrepr(self):
        """Convert to V-representation.

        This potentially expensive operation is implemented by calling the ``lrs`` tool from ``lrslib``.

        :rtype: :class:`VRepr`"""
        stdout = _call_lrs('lrs', self._to_lrs())
        return VRepr._from_lrs(stdout)

    def irred(self):
        """Return H-representation with redundant inequalities and equations removed.

        *Warning:* This is unique only if there are no equations.

        This potentially expensive operation is implemented by calling the ``redund`` tool from ``lrslib``.

        :rtype: :class:`HRepr`"""
        stdout = _call_lrs('redund', self._to_lrs())
        return HRepr._from_lrs(stdout)

    def vertices(self):
        """Return vertices of convex polytope described by this H-representation.

        Convenience function that ensures that polyhedron is bounded.
        This potentially expensive operation is implemented by calling :func:`vrepr`.

        :rtype: list of :class:`vector`
        :raises ValueError: polytope is not bounded (use :func:`vrepr` to access V-representation)
        """
        vrepr = self.vrepr()
        if len(vrepr.rays) > 0 or len(vrepr.lines) > 0:
            raise ValueError, 'Polyhedron has extremal rays.'
        return vrepr.vertices

    def __and__(self, rhs):
        """Intersect two H-representations (the resulting H-representation will typically be redundant).

        :param rhs: the H-representation to intersect with.
        :rtype: :class:`HRepr`
        """
        assert self.ambient_dim == rhs.ambient_dim
        return HRepr(
            ieqs=self.ieqs + rhs.ieqs,
            eqns=self.eqns + rhs.eqns,
            ambient_dim=self.ambient_dim)

    # def map(self, f):
    #     return HRepr(
    #         ieqs=map(f, self.ieqs),
    #         eqns=map(f, self.eqns),
    #         ambient_dim=self.ambient_dim)

    def to_sage(self):
        """Convert to Sage :class:`Polyhedron` object.

        :rtype: :class:`sage.Polyhedron`"""
        if self.eqns or self.ieqs:
            logger.debug('Converting HRepr to Sage Polyhedron')

            def convert((H, c)):
                return [-c] + list(H)

            ieqs = map(convert, self.ieqs)
            eqns = map(convert, self.eqns)

            return Polyhedron(ieqs=ieqs, eqns=eqns)
        else:
            return Polyhedron(eqns=[[0] * (1 + self.ambient_dim)])

    def _to_lrs(self):
        """Convert H-representation to ``lrs`` format.

        :rtype: string"""
        lines = ['my_polytope']
        lines += ['H-representation']
        if self.eqns:
            lines += [
                'linearity %d %s' %
                (len(self.eqns),
                 ' '.join(map(str, range(1, len(self.eqns) + 1))))
            ]
        lines += ['begin']
        lines += [
            '%d %d rational' %
            (len(self.eqns) + len(self.ieqs), 1 + self.ambient_dim)
        ]
        lines += [
            str(-c) + ' ' + ' '.join(map(str, H)) for (H, c) in self.eqns
        ]
        lines += [
            str(-c) + ' ' + ' '.join(map(str, H)) for (H, c) in self.ieqs
        ]
        lines.append('end')
        return '\n'.join(lines)

    @staticmethod
    def from_sage(p):
        """Construct H-representation from Sage :class:`Polyhedron`.

        :param p: the :class:`Polyhedron` instance.
        :rtype: :class:`HRepr`
        """

        def convert(ieq):
            return (ieq[1:], -ieq[0])

        ieqs = map(convert, p.inequalities())
        eqns = map(convert, p.equations())
        return HRepr(ieqs=ieqs, eqns=eqns)

    @staticmethod
    def _from_lrs(s):
        """Construct H-representation from `lrs` output."""
        # split into lines and ignore comments
        lines = [
            line.strip() for line in s.splitlines()
            if line.strip() and (line.startswith('*****') or
                                 not line.startswith('*'))
        ]

        # parse header
        name = lines.pop(0)
        type = lines.pop(0)
        assert name == 'my_polytope'
        assert type == 'H-representation'

        # linearity?
        if lines[0].startswith('linearity'):
            lin = lines.pop(0).split()
            assert lin[0] == 'linearity'
            lin = map(int, lin[1:])
            assert len(lin) == lin[0] + 1
            lin = lin[1:]
        else:
            lin = []

        # (number of inequalities or *****) and dimension
        assert lines.pop(0) == 'begin'
        rat = lines.pop(0).split()
        assert rat[2] == 'rational'
        #rat = map(int, rat[:2])

        ieqs = []
        eqns = []
        for i, line in enumerate(lines):
            if line == 'end': break
            line = map(QQ, line.split())
            c = -line[0]
            H = vector(line[1:])
            if (i + 1) in lin:
                eqns.append((H, c))
            else:
                ieqs.append((H, c))

        return HRepr(ieqs=ieqs, eqns=eqns)


class VRepr(object):
    """Light-weight container that stores a polyhedron in V-representation.

    A V-representation without any vertices, rays, and lines is assumed to be empty.
    Vertices, rays, and lines are sorted lexicographically.

    :param vertices: vertices.
    :param rays: rays.
    :param lines: lines.
    :type vertices: list of :class:`vector`
    :type rays: list of :class:`vector`
    :type lines: list of :class:`vector`
    :raises ValueError: neither ``vertices``, ``rays``, ``lines``, nor ``ambient_dim`` were specified.

    .. automethod:: __eq__
    """

    def __init__(self, vertices=[], rays=[], lines=[], ambient_dim=None):
        # sort vertices, rays, lines lexicographically
        #: The vertices.
        self.vertices = sorted(map(vector, vertices))
        #: The rays.
        self.rays = sorted(map(vector, rays))
        #: The lines.
        self.lines = sorted(map(vector, lines))

        # determine ambient dimension
        if not (ambient_dim or (vertices or rays or lines)):
            raise ValueError(
                'Need to specify ambient dimension if vertices, rays, and lines are empty.'
            )
        if not ambient_dim:
            V = next(
                iter(
                    self.vertices if self.vertices else (self.rays if self.rays
                                                         else self.lines)))
            ambient_dim = len(V)
        #: The ambient dimension.
        self.ambient_dim = ambient_dim

        # check format
        assert all(len(V) == ambient_dim for V in self.vertices)
        assert all(len(V) == ambient_dim for V in self.rays)
        assert all(len(V) == ambient_dim for V in self.lines)

    def __str__(self):
        return '<VRepr object at 0x%016X in dimension %d defined by %d vertices, %d rays and %d lines>' % (
            id(self), self.ambient_dim, len(self.vertices), len(self.rays),
            len(self.lines))

    def __eq__(self, rhs):
        """Compare V-representations.

        *Warning:* Equal polyhedra can have different V-representation!

        :rtype: bool"""
        return self.vertices == rhs.vertices and self.rays == rhs.rays and self.lines == rhs.lines

    def hrepr(self):
        """Convert to H-representation.

        This potentially expensive operation is implemented by calling the ``lrs`` tool from ``lrslib``.

        :rtype: :class:`HRepr`"""
        stdout = _call_lrs('lrs', self._to_lrs())
        return HRepr._from_lrs(stdout)

    def to_sage(self):
        """Convert to Sage :class:`Polyhedron` object.

        :rtype: :class:`sage.Polyhedron`"""
        logger.debug('Converting VRepr to Sage Polyhedron')

        return Polyhedron(
            vertices=self.vertices,
            rays=self.rays,
            lines=self.lines,
            ambient_dim=self.ambient_dim)

    def _to_lrs(self):
        """Convert V-representation to ``lrs`` format.

        :rtype: string"""
        lines = ['my_polytope']
        lines += ['V-representation']
        if self.lines:
            lines += [
                'linearity %d %s' %
                (len(self.lines),
                 ' '.join(map(str, range(1, len(self.lines) + 1))))
            ]
        lines += ['begin']
        lines += [
            '%d %d rational' %
            (len(self.lines) + len(self.rays) + len(self.vertices),
             1 + self.ambient_dim)
        ]
        lines += ['0 ' + ' '.join(map(str, V)) for V in self.lines]
        lines += ['0 ' + ' '.join(map(str, V)) for V in self.rays]
        lines += ['1 ' + ' '.join(map(str, V)) for V in self.vertices]
        lines.append('end')
        return '\n'.join(lines)

    @staticmethod
    def from_sage(p):
        """Construct V-representation from Sage :class:`Polyhedron`.

        :param p: the :class:`Polyhedron` instance.
        :rtype: :class:`VRepr`
        """
        return VRepr(vertices=p.vertices(), rays=p.rays(), lines=p.lines())

    @staticmethod
    def _from_lrs(s):
        """Construct V-representation from `lrs` output."""
        lines = [
            line.strip() for line in s.splitlines()
            if line.strip() and (line.startswith('*****') or
                                 not line.startswith('*'))
        ]

        # parse header
        name = lines.pop(0)
        type = lines.pop(0)
        assert name == 'my_polytope'
        assert type == 'V-representation'

        # extremal rays?
        if lines[0].startswith('linearity'):
            lin = lines.pop(0).split()
            assert lin[0] == 'linearity'
            lin = map(int, lin[1:])
            assert len(lin) == lin[0] + 1
            lin = lin[1:]
        else:
            lin = []

        # (number of entries or *****) and dimension
        assert lines.pop(0) == 'begin'
        rat = lines.pop(0).split()
        assert rat[2] == 'rational'
        #rat = map(int, rat[:2])

        lins = []
        rays = []
        vertices = []
        for i, line in enumerate(lines):
            if line == 'end': break
            line = line.split()

            # parse line
            type = int(line[0])
            line = vector(map(QQ, line[1:]))

            # vertex or extremal ray?
            if type == 0:
                if (i + 1) in lin:
                    lins.append(line)
                else:
                    rays.append(line)
            elif type == 1:
                vertices.append(line)
            else:
                assert type in (0, 1), 'Unknown type: %d' % type
        return VRepr(vertices=vertices, rays=rays, lines=lins)
