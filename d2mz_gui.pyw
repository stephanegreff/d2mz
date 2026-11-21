#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""d2mz_gui — point-and-click front-end for d2mz (standard library only).

Double-click launcher (see d2mz.bat on Windows, or run this file with
pythonw): add one or several folders containing your Bruker .d
acquisitions to the conversion queue (Add folder... / drag & drop /
type a path + Enter / Remove selected), optionally check "Single output folder" and pick a
destination shared by the whole queue, pick the output format, press
Convert. The whole queue is converted in order. By default each
source folder gets its own d2mz_output_<source-name>/ folder created
NEXT TO it, with the mzML and mzXML files in its analyses/ subfolder
and the precursor tables in its annexes/TSV/ subfolder, while the
series report_<source>.json and the CITATION files sit in its
annexes/ subfolder. With "Single output folder" everything goes to the
chosen path instead.

No dependency beyond the Python standard library (tkinter included).

2026-11-09 (v_bi): the conversion worker now logs through the
main-thread queue (tkinter is not thread-safe; direct Text updates
from the worker crash silently under pythonw and made Convert look
like a dead button), the Convert handler traces its entry and logs
full tracebacks into the log window, and argparse rejections
(SystemExit) are reported instead of killing the worker silently.

2026-11-09 (v_bj): WM_DROPFILES fixed to 0x0233 (the previous
0x0232 never matched, so drops were silently ignored), and the
bottom buttons reordered on S. request: Open on the left, Quit in
the middle, Convert on the right.
"""
import os
import sys
import threading
import queue
import contextlib
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d2mz

try:
    import tkinter as tk
    from tkinter import filedialog, ttk
except ImportError:            # headless machine: allow import for testing
    tk = None


def build_argvec(folder, fmt, precision, out=None, gz=False,
                  zipout=False):
    """Command-line equivalent of the GUI settings (used by tests).
    Polarity is detected automatically by d2mz itself. `out` is the
    optional single output folder for the whole queue (--out);
    `zipout` adds --zip-output (archive the output folder)."""
    argvec = [folder, '--format', fmt, '--precision', precision]
    if gz:
        argvec.append('--gzip')
    if zipout:
        argvec.append('--zip-output')
    if out:
        argvec += ['--out', out]
    return argvec


# ------------------------------------------------- Windows drag & drop
def _enable_drop(root, add_paths, log=None):
    """Drag & drop of folders/.zip archives onto the d2mz window on
    Windows, with the standard library only (ctypes WM_DROPFILES
    hook, the windnd technique). No-op on other platforms."""
    import ctypes
    import ctypes.wintypes as wintypes
    import struct as _struct
    user32 = ctypes.windll.user32
    shell32 = ctypes.windll.shell32
    GWL_WNDPROC = -4
    # v_bj fix (drop never received): WM_DROPFILES is 0x0233 in
    # WinUser.h. The previous 0x0232 never matched, so every drop was
    # silently forwarded to the default window procedure and ignored
    WM_DROPFILES = 0x0233
    WM_DESTROY = 0x0002
    # WNDPROC does NOT exist in the ctypes module: the callback
    # prototype is defined here (windnd technique). LRESULT is a
    # pointer-sized signed integer.
    if _struct.calcsize('P') * 8 == 64:
        _LRESULT = ctypes.c_longlong
    else:
        _LRESULT = ctypes.c_long
    WNDPROC = ctypes.WINFUNCTYPE(
        _LRESULT, wintypes.HWND, ctypes.c_uint,
        wintypes.WPARAM, wintypes.LPARAM)
    # 32-bit Windows has no SetWindowLongPtrW/GetWindowLongPtrW:
    # the Long aliases do the same job there.
    try:
        _set_wndproc = user32.SetWindowLongPtrW
        _get_wndproc = user32.GetWindowLongPtrW
    except AttributeError:
        _set_wndproc = user32.SetWindowLongW
        _get_wndproc = user32.GetWindowLongW
    # fix 2 (crash 88): WITHOUT argtypes the 64-bit callback pointer
    # overflows the default C int argument. WNDPROC pointers travel
    # as c_void_p (the function-pointer instance is accepted as is).
    _get_wndproc.argtypes = [wintypes.HWND, ctypes.c_int]
    _get_wndproc.restype = ctypes.c_void_p
    _set_wndproc.argtypes = [wintypes.HWND, ctypes.c_int, ctypes.c_void_p]
    _set_wndproc.restype = ctypes.c_void_p
    user32.CallWindowProcW.argtypes = [
        ctypes.c_void_p, wintypes.HWND, ctypes.c_uint,
        wintypes.WPARAM, wintypes.LPARAM]
    user32.CallWindowProcW.restype = _LRESULT
    # fix 3 (crash 89): DragQueryFileW/DragFinish live in shell32,
    # NOT in user32 (windnd uses shell32 for both).
    shell32.DragQueryFileW.argtypes = [
        wintypes.HDROP, ctypes.c_uint, ctypes.c_wchar_p, ctypes.c_uint]
    shell32.DragQueryFileW.restype = ctypes.c_uint
    shell32.DragFinish.argtypes = [wintypes.HDROP]
    # DragAcceptFiles is exported by shell32 as well (the whole
    # WM_DROPFILES family lives there)
    shell32.DragAcceptFiles.argtypes = [wintypes.HWND, wintypes.BOOL]
    # fix 4 (crash report 89, drop never received): register the
    # drop on the REAL top-level window, however deep winfo_id sits
    # in the Tk wrapper chain (walk the whole parent chain up).
    # v_cd (drop refused -- forbidden cursor): make_gui installs the
    # hook BEFORE the Tk window is realized, so winfo_id() can return
    # a not-yet-created/invalid native window and DragAcceptFiles then
    # registers nothing: the cursor stays "forbidden". Realize the
    # window first, and register on BOTH the root native window and
    # the real top-level ancestor (whichever the mouse hovers, the
    # drop lands on one of them).
    root.update_idletasks()
    root.update()
    _child = wintypes.HWND(root.winfo_id())
    hwnd = _child
    _n_up = 0
    while True:
        parent = user32.GetParent(hwnd)
        if not parent:
            break
        hwnd = parent
        _n_up += 1
    if log:
        log('drop hook: root native hwnd=%s, top-level hwnd=%s '
            '(%d level(s) up)' % (hex(_child.value or 0),
                                  hex(hwnd.value or 0), _n_up))
        if not (hwnd.value or 0):
            log('drop hook ERROR: null hwnd -- window not realized; '
                'drop registration aborted')

    class _Handler:
        def __init__(self):
            # v_cd: swap the wndproc on the TOP-LEVEL window (it gets
            # every WM_DROPFILES), but ALSO mark the root native
            # window as accepting files -- on some Tk builds the
            # visible surface is the child, and Windows shows the
            # "forbidden" cursor when the hovered window itself does
            # not accept drops.
            self.old = _get_wndproc(hwnd, GWL_WNDPROC)
            self.proc = WNDPROC(self._wndproc)
            _set_wndproc(hwnd, GWL_WNDPROC,
                         ctypes.cast(self.proc, ctypes.c_void_p))
            shell32.DragAcceptFiles(hwnd, True)
            if _child.value and _child.value != hwnd.value:
                shell32.DragAcceptFiles(_child, True)

        def _wndproc(self, h, msg, wparam, lparam):
            if msg == WM_DROPFILES:
                try:
                    hdrop = wintypes.HDROP(wparam)
                    n = shell32.DragQueryFileW(hdrop, 0xFFFFFFFF,
                                               None, 0)
                    buf = ctypes.create_unicode_buffer(1024)
                    paths = []
                    for i in range(n):
                        shell32.DragQueryFileW(hdrop, i, buf, 1024)
                        paths.append(buf.value)
                    shell32.DragFinish(hdrop)
                    if log:
                        log('drop received: ' + ', '.join(paths))
                    root.after(0, lambda: add_paths(paths))
                except Exception as e:
                    if log:
                        log('drop handling failed: %r' % (e,))
                return 0
            if msg == WM_DESTROY:
                _set_wndproc(hwnd, GWL_WNDPROC,
                             ctypes.c_void_p(self.old or 0))
            return user32.CallWindowProcW(self.old, h, msg,
                                          wparam, lparam)

    return _Handler()


# --------------------------------------------------------------------- GUI
def make_gui(root):
    state = {'worker': None, 'queue': queue.Queue()}

    root.title('d2mz — Bruker .d → mzML/mzXML')
    root.minsize(640, 420)
    root.geometry('780x600')
    frm = ttk.Frame(root, padding=10)
    frm.pack(fill='both', expand=True)

    # --- conversion queue: ONE OR SEVERAL folders (v_ay: listbox +
    # Add/Remove; the Convert pass processes the queue in order)
    ttk.Label(frm, text='Folder(s) or .zip archive(s) with .d '
              'acquisitions (queue — several entries allowed):').grid(
        row=0, column=0, columnspan=3, sticky='w')
    var_folder = tk.StringVar()
    ent = ttk.Entry(frm, textvariable=var_folder, width=60)
    ent.grid(row=1, column=0, columnspan=2, sticky='we', padx=(0, 5))

    queue_box = tk.Listbox(frm, height=4, selectmode='extended',
                           activestyle='dotbox')
    queue_box.grid(row=2, column=0, columnspan=2, sticky='nsew',
                   pady=(4, 0))

    def add_folder(d):
        d = (d or '').strip()
        ok = (os.path.isdir(d) or
              (os.path.isfile(d) and d.lower().endswith('.zip')))
        if ok and d not in queue_box.get(0, 'end'):
            queue_box.insert('end', d)
    def browse_folder():
        # analyst request 2026-11-07: the dialog REOPENS after each pick
        # (Cancel to stop) so that several folders are queued in a row
        while True:
            d = filedialog.askdirectory(
                title='Folder with .d acquisitions (Cancel when done)')
            if not d:
                break
            add_folder(d)
        var_folder.set('')
    def browse_zip():
        # v_az: .zip archives of .d acquisitions are accepted too
        d = filedialog.askopenfilename(
            title='.zip archive of .d acquisitions',
            filetypes=[('Zip archive', '*.zip'), ('All files', '*.*')])
        if d:
            add_folder(d)
            var_folder.set('')
    def run_probe():
        # v_ce (analyst request 2026-11-12 -- zero console): diagnose a
        # .d with d2mz's probe and write the report NEXT TO the .d,
        # then open that folder so the report is one drag away.
        d = filedialog.askdirectory(
            title='Diagnostic: select ONE .d folder (e.g. Isolated_A1...37320.d)')
        if not d:
            return
        d = d.rstrip('\\/')
        if not d.lower().endswith('.d'):
            log_line('diagnostic: %s is not a .d folder -- pick one '
                     'of the folders ending with .d' % d)
            return
        log_line('')
        log_line('diagnostic of %s (this can take ~1 min) ...' % d)
        lines = []
        sweep_lines = []

        def plog(msg):
            lines.append(str(msg))
            log_line(str(msg))
        try:
            d2mz.probe_folder(d, plog, sweep_log=sweep_lines.append)
            rpt = os.path.join(
                os.path.dirname(d), 'd2mz_probe_%s_%s.txt'
                % (os.path.basename(d)[:-2], d2mz.__version__))
            with open(rpt, 'w', encoding='utf-8') as f:
                f.write('\n'.join(lines) + '\n')
            log_line('diagnostic done -> report saved:')
            log_line('    ' + rpt)
            # v_cf: the SWEEP dump goes to a second, SMALL file -- the
            # full report exceeds the 4 MiB chat upload limit
            if sweep_lines:
                rpt2 = os.path.join(
                    os.path.dirname(d), 'd2mz_probe_sweep_%s_%s.txt'
                    % (os.path.basename(d)[:-2], d2mz.__version__))
                with open(rpt2, 'w', encoding='utf-8') as f:
                    f.write('\n'.join(sweep_lines) + '\n')
                log_line('sweep-MS2 report saved (small, uploadable):')
                log_line('    ' + rpt2)
                log_line('Please upload the SWEEP report '
                         '(d2mz_probe_sweep_...txt).')
            try:
                os.startfile(os.path.dirname(d))
            except Exception:
                pass
        except Exception:
            log_line('diagnostic FAILED -- please report this block:')
            for line in traceback.format_exc().splitlines():
                log_line(line)

    def add_typed(_evt=None):
        add_folder(var_folder.get())
        var_folder.set('')
    def remove_selected():
        for i in sorted(queue_box.curselection(), reverse=True):
            queue_box.delete(i)

    btns = ttk.Frame(frm)
    btns.grid(row=1, column=2, rowspan=2, sticky='nw', padx=(5, 0))
    ttk.Button(btns, text='Add folder…',
               command=browse_folder).pack(fill='x', pady=(0, 3))
    ttk.Button(btns, text='Add .zip…',
               command=browse_zip).pack(fill='x', pady=(0, 3))
    ttk.Button(btns, text='Diagnostic (probe)…',
               command=run_probe).pack(fill='x', pady=(0, 3))
    ttk.Button(btns, text='Remove selected',
               command=remove_selected).pack(fill='x')
    ent.bind('<Return>', add_typed)

    # --- single output folder for the whole queue (optional; default
    # since v_bd: each source folder gets its own d2mz_output_<name>/
    # folder created NEXT TO it, so several sources never mix)
    var_single_out = tk.BooleanVar(value=False)
    var_out = tk.StringVar()
    ttk.Checkbutton(frm, text='Single output folder for the whole queue',
                   variable=var_single_out).grid(
        row=3, column=0, columnspan=3, sticky='w', pady=(10, 0))
    ent_out = ttk.Entry(frm, textvariable=var_out, width=60)
    ent_out.grid(row=4, column=0, columnspan=2, sticky='we',
                padx=(16, 5), pady=(2, 0))

    def browse_out():
        d = filedialog.askdirectory(title='Single output folder')
        if d:
            var_out.set(d)
            var_single_out.set(True)

    def _sync_out_state():
        normal = 'normal' if var_single_out.get() else 'disabled'
        ent_out.configure(state=normal)
    var_single_out.trace_add('write', lambda *_: _sync_out_state())

    ttk.Button(frm, text='Browse…', command=browse_out).grid(
        row=4, column=2, sticky='w', pady=(2, 0))
    ent_out.configure(state='disabled')
    ent_out.bind('<Return>', lambda _e: var_single_out.set(True))

    ttk.Label(frm, text='Output format:').grid(
        row=5, column=0, columnspan=3, sticky='w', pady=(10, 0))
    var_fmt = tk.StringVar(value='both')
    for i, (lbl, val) in enumerate((('mzML + mzXML', 'both'),
                                    ('mzML only', 'mzml'),
                                    ('mzXML only', 'mzxml'))):
        ttk.Radiobutton(frm, text=lbl, variable=var_fmt,
                        value=val).grid(row=6, column=i, sticky='w')

    var_prec = tk.StringVar(value='32')
    ttk.Label(frm, text='Binary precision (bits):').grid(
        row=7, column=0, sticky='w', pady=(10, 0))
    for i, b in enumerate(('32', '64')):
        ttk.Radiobutton(frm, text=b, variable=var_prec,
                        value=b).grid(row=7, column=1 + i, sticky='w')

    var_gz = tk.BooleanVar(value=False)
    var_zipout = tk.BooleanVar(value=False)
    ttk.Checkbutton(frm, variable=var_gz,
                    text='Gzip converted files (.mzML.gz / .mzXML.gz, '
                         'lossless, read directly by MZmine)'
                    ).grid(row=8, column=0, columnspan=3, sticky='w',
                           pady=(6, 0))
    # v_bl: transfer/archiving option. MZmine still reads the flat
    # .mzML/.mzML.gz files; a .zip must be extracted first, which is
    # why the flat files are always kept.
    ttk.Checkbutton(frm, variable=var_zipout,
                    text='Also write the output folder as ONE .zip '
                         'archive next to it (transfer/archiving; '
                         'the flat files stay for MZmine)'
                    ).grid(row=9, column=0, columnspan=3, sticky='w',
                           pady=(6, 0))

    var_status = tk.StringVar(
        value='Polarity is detected automatically for every run; '
              'lock-mass recalibration uses the built-in Bruker sodium '
              'formate/acetate cluster series.')
    ttk.Label(frm, textvariable=var_status, wraplength=600,
              foreground='#555').grid(row=10, column=0, columnspan=3,
                                      sticky='w', pady=(8, 2))

    txt = tk.Text(frm, height=10, width=80, state='disabled')
    txt.grid(row=11, column=0, columnspan=3, sticky='nsew', pady=5)

    def _log_insert(line):
        # touches the Text widget — main thread ONLY (v_bi:
        # tkinter calls from the worker thread crash silently
        # under pythonw, which made Convert look like a dead button)
        txt.configure(state='normal')
        txt.insert('end', line + '\n')
        txt.see('end')
        txt.configure(state='disabled')

    def log_line(line):
        # thread-safe log: the conversion worker never touches
        # tkinter directly, it hands its lines to the queue and the
        # main thread drains it in poll()
        if threading.current_thread() is threading.main_thread():
            _log_insert(line)
        else:
            state['queue'].put(line)

    def poll():
        while True:
            try:
                item = state['queue'].get_nowait()
            except queue.Empty:
                break
            if item is True:            # end-of-run marker (kept for
                _log_insert('')         # older d2mz_gui copies)
            else:
                _log_insert(item)
        root.after(150, poll)

    root.after(150, poll)

    def run_conversion():
        # v_bi: diagnose "Convert does nothing" reports. The button
        # command MUST leave a trace in the log whatever happens next:
        # pythonw swallows stderr, so an exception in this handler
        # would otherwise look like a dead button.
        log_line('Convert: clicked.')
        try:
            _run_conversion_inner()
        except Exception:
            log_line('ERROR in Convert handler — please report '
                     'this block:')
            for line in traceback.format_exc().splitlines():
                log_line(line)

    def _run_conversion_inner():
        folders = list(queue_box.get(0, 'end'))
        if not folders:
            log_line('!! Add at least one folder to the queue first.')
            return
        bad = [f for f in folders
               if not (os.path.isdir(f)
                       or (os.path.isfile(f)
                           and f.lower().endswith('.zip')))]
        if bad:
            for f in bad:
                log_line('!! Not a folder: %s' % f)
            return
        if state['worker'] is not None and state['worker'].is_alive():
            return
        fmt, prec = var_fmt.get(), var_prec.get()
        out = None
        if var_single_out.get():
            out = var_out.get().strip()
            if not out:
                log_line('!! Check "Single output folder" needs a path.')
                return
            try:
                os.makedirs(out, exist_ok=True)
            except OSError as e:
                log_line('!! Cannot create output folder: %s' % e)
                return

        class _Sink:
            def write(self, s):
                if s.strip():
                    for part in s.splitlines():
                        log_line(part.rstrip())
            def flush(self):
                pass

        def work():
            rcs = []   # v_cb: per-folder d2mz exit codes, surfaced
                       # at the end so an error can never hide
                       # behind a bare "process finished."
            tot_warn = 0   # v_cc (analyst request 2026-11-12): quality
                           # warnings (e.g. 0-MS2 baf-fallback files)
                           # must end on an explicit line too, not on
                           # a bare "process finished. no error."
            try:
                with contextlib.redirect_stdout(_Sink()):
                    for folder in folders:
                        argvec = build_argvec(folder, fmt, prec, out,
                                               var_gz.get(),
                                               var_zipout.get())
                        log_line('')
                        log_line('$ python3 d2mz.py ' + ' '.join(argvec))
                        try:
                            rcs.append(d2mz.main(argvec))
                            tot_warn += getattr(
                                d2mz, 'LAST_RUN',
                                {}).get('n_warnings', 0)
                        except SystemExit as e:  # argparse rejection
                            rcs.append(e.code)
                            log_line('ERROR: d2mz exited with code %s '
                                     '(bad argument?)' % e.code)
                        except Exception:        # pragma: no cover
                            rcs.append(1)
                            log_line('ERROR inside d2mz — please report '
                                     'this block:')
                            for line in (traceback.format_exc()
                                         ).splitlines():
                                log_line(line)
            except Exception:                   # pragma: no cover
                rcs.append(1)
                log_line('ERROR in conversion worker — please report '
                         'this block:')
                for line in traceback.format_exc().splitlines():
                    log_line(line)
            finally:
                if any(rc != 0 and rc is not None for rc in rcs):
                    log_line('process finished WITH ERRORS '
                             '(exit codes: %s) -- output '
                             'INCOMPLETE, see the recap above.'
                             % ', '.join(str(rc) for rc in rcs))
                elif tot_warn:
                    log_line('process finished -- %d QUALITY '
                             'WARNING(S): conversion DEGRADED (e.g. '
                             'files without MS2), see the QUALITY '
                             'WARNINGS list above.' % tot_warn)
                else:
                    log_line('process finished. no error.')

        state['worker'] = threading.Thread(target=work)
        state['worker'].daemon = True
        state['worker'].start()

    def open_output():
        if var_single_out.get() and var_out.get().strip():
            out = var_out.get().strip()
        else:
            sel = queue_box.curselection()
            folder = (queue_box.get(sel[0]) if sel else
                      (queue_box.get(0) if queue_box.size() else ''))
            # v_bd: default output sits NEXT TO the source folder,
            # named after it (d2mz_output_<source-name>/)
            if folder:
                src = os.path.normpath(folder)
                out = os.path.join(os.path.dirname(src),
                                  'd2mz_output_' + os.path.basename(src))
            else:
                out = ''
        if out and os.path.isdir(out):
            if sys.platform == 'win32':
                os.startfile(out)              # noqa: attribution windows
            elif sys.platform == 'darwin':
                os.system('open "%s"' % out)
            else:
                os.system('xdg-open "%s"' % out)

    # 2026-11-09 (S. request): Open on the LEFT (he kept clicking it
    # instead of Convert), Quit in the middle, Convert on the RIGHT
    ttk.Button(frm, text='Open d2mz_output folder',
               command=open_output).grid(row=12, column=0, sticky='w',
                                         pady=5)
    ttk.Button(frm, text='Quit', command=root.destroy).grid(
        row=12, column=1, pady=5)
    ttk.Button(frm, text='Convert', command=run_conversion).grid(
        row=12, column=2, sticky='e', pady=5)
    frm.columnconfigure(0, weight=1)
    # rows 0-2 (queue) and 10 (log) stretch with the window; the
    # settings block in between keeps its fixed height
    frm.rowconfigure(2, weight=1)
    frm.rowconfigure(11, weight=2)
    # drag & drop AFTER log_line exists: the journal must say whether
    # the drop hook actually got installed and what it receives
    if sys.platform == 'win32':
        def _install_drop():
            # v_cd: install only once the window is REALIZED and
            # mapped (deferred to after the first layout pass -- the
            # naive install at make_gui time could read winfo_id()
            # before the native window existed and register drops on
            # a null handle: forbidden cursor, drop never received).
            try:
                # the handler MUST stay referenced: losing the WNDPROC
                # callback to the garbage collector would crash the
                # window
                state['drop'] = _enable_drop(
                    root, lambda paths: [add_folder(p) for p in paths],
                    log_line)
                log_line('drag & drop enabled: folders and .zip '
                         'files can be dropped on this window')
            except Exception as e:
                log_line('drag & drop unavailable (%r) — use Add '
                         'folder or Add .zip instead' % (e,))
        root.update_idletasks()
        root.after(150, _install_drop)
    return state


def main():
    root = tk.Tk()
    make_gui(root)
    root.mainloop()


if __name__ == '__main__':
    if tk is None:                                  # headless fallback
        print('tkinter unavailable — use: python3 d2mz.py <folder>')
        sys.exit(2)
    try:
        main()
    except Exception:
        import traceback
        tb = traceback.format_exc()
        try:                                        # keep a log next to the script
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'd2mz_gui_error.log'), 'w') as f:
                f.write(tb)
        except OSError:
            pass
        try:
            from tkinter import messagebox
            messagebox.showerror('d2mz — startup error', tb)
        except Exception:
            pass
        sys.exit(1)
