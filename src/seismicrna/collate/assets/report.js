/* SEISMIC-RNA collated report — dashboard behavior.
 *
 * Depends on globals injected by the template:
 *   window.__reportConfig  = { reportId, groupIds: [...] }
 *   window.__seismicPlotly = { "<version>": PlotlyModule, ... }
 *   window.__seismicRender = [ function(){...}, ... ]  (plot render thunks)
 *   GridStack (from gridstack-all.js)
 *   svgPanZoom (optional, from svg-pan-zoom.js)
 */
(function () {
    "use strict";

    var CFG = window.__reportConfig || { reportId: "report", groupIds: [] };
    // Bump this version to invalidate previously saved layouts when the
    // default layout scheme changes.
    var STORE_PREFIX = "seismic-collate/v5/";
    var LAYOUT_KEY = function (gid) { return STORE_PREFIX + CFG.reportId + "/layout/" + gid; };
    var COLLAPSE_KEY = function (gid) { return STORE_PREFIX + CFG.reportId + "/collapsed/" + gid; };
    var ORDER_KEY = STORE_PREFIX + CFG.reportId + "/grouporder";
    var THEME_KEY = STORE_PREFIX + "theme";

    var grids = {};          // groupId -> GridStack instance
    var defaultLayouts = {}; // groupId -> layout snapshot (positions only)
    var defaultGroupOrder = null; // group ids in their original order
    var suppressPersist = false; // true while applying layout programmatically
    var panZoomInstances = [];

    // Run a layout change that must NOT be written back to localStorage
    // (e.g. loading a saved layout on open).
    function programmatic(fn) {
        suppressPersist = true;
        try { fn(); } finally {
            setTimeout(function () { suppressPersist = false; }, 450);
        }
    }

    function store(key, value) {
        try { localStorage.setItem(key, JSON.stringify(value)); } catch (e) { /* quota / privacy mode */ }
    }
    function load(key) {
        try {
            var raw = localStorage.getItem(key);
            return raw ? JSON.parse(raw) : null;
        } catch (e) { return null; }
    }
    function remove(key) {
        try { localStorage.removeItem(key); } catch (e) { /* ignore */ }
    }

    /* ---------------- Theme ---------------- */
    function applyTheme(theme) {
        var root = document.documentElement;
        if (theme === "light" || theme === "dark") {
            root.setAttribute("data-theme", theme);
        } else {
            // auto: follow system
            var dark = window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
            root.setAttribute("data-theme", dark ? "dark" : "light");
        }
        var btn = document.getElementById("theme-btn");
        if (btn) {
            var eff = root.getAttribute("data-theme");
            btn.querySelector(".label").textContent = theme === "auto" ? "Auto" : (eff === "dark" ? "Dark" : "Light");
        }
    }
    function initTheme() {
        var theme = load(THEME_KEY) || "auto";
        applyTheme(theme);
        var btn = document.getElementById("theme-btn");
        if (btn) {
            btn.addEventListener("click", function () {
                var cur = load(THEME_KEY) || "auto";
                var next = cur === "auto" ? "light" : (cur === "light" ? "dark" : "auto");
                store(THEME_KEY, next);
                applyTheme(next);
                resizeAllPlots();
            });
        }
        if (window.matchMedia) {
            window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", function () {
                if ((load(THEME_KEY) || "auto") === "auto") { applyTheme("auto"); resizeAllPlots(); }
            });
        }
    }

    /* ---------------- Lazy tile activation ----------------
     * Rendering every plot (portable) or loading every graph iframe (linked)
     * up front is wasteful when most tiles are off-screen. Each tile is
     * activated only once it nears the viewport: the inline plot is rendered,
     * or the linked iframe's src is set. */
    function activateTile(item) {
        if (!item || item._activated) return;
        item._activated = true;
        // Inline (portable) Plotly plot: run its render thunk once.
        var div = item.querySelector(".plotly-graph-div");
        var map = window.__seismicRender;
        if (div && map && map[div.id]) {
            var fn = map[div.id];
            delete map[div.id];
            try { fn(); } catch (e) {
                console.error("SEISMIC report: failed to render a plot", e);
            }
        }
        // Linked (non-portable) graph iframe: load it now.
        var frame = item.querySelector("iframe[data-src]");
        if (frame && !frame.src) {
            frame.addEventListener("load", function () {
                setTimeout(function () { resizePlotsIn(item); }, 30);
            });
            frame.src = frame.getAttribute("data-src");
        }
    }
    function initLazyTiles() {
        var items = document.querySelectorAll(".grid-stack-item");
        if (typeof IntersectionObserver === "undefined") {
            items.forEach(activateTile);
            return;
        }
        // Activate a bit before a tile scrolls into view for a seamless feel.
        var io = new IntersectionObserver(function (entries) {
            entries.forEach(function (en) {
                if (en.isIntersecting) {
                    activateTile(en.target);
                    io.unobserve(en.target);
                }
            });
        }, { root: null, rootMargin: "600px 0px", threshold: 0 });
        items.forEach(function (it) { io.observe(it); });
        window.__tileIO = io;
    }

    // Resize the inline/iframe Plotly plots in `el` to fill their tiles. This is
    // safe to call from a ResizeObserver on every tile size change; it never
    // touches pan/zoom drawings (whose resize mode depends on *why* the size
    // changed, so those are driven only by the explicit handlers below).
    function resizeFramedPlots(el) {
        if (!el) return;
        el.querySelectorAll(".plotly-graph-div").forEach(function (d) {
            // Only rendered plots have _seismicPlotly; skip ones not yet activated.
            var mod = d._seismicPlotly;
            if (mod && mod.Plots && d.offsetParent !== null) {
                try { mod.Plots.resize(d); } catch (e) { /* ignore */ }
            }
        });
        el.querySelectorAll("iframe[data-resizable]").forEach(function (f) {
            try {
                var w = f.contentWindow;
                if (!w) return;
                var plots = w.document.querySelectorAll(".plotly-graph-div");
                if (w.Plotly && w.Plotly.Plots && plots.length) {
                    plots.forEach(function (d) {
                        try { w.Plotly.Plots.resize(d); } catch (e) { /* ignore */ }
                    });
                } else {
                    w.dispatchEvent(new Event("resize"));
                }
            } catch (e) { /* iframe not yet loaded */ }
        });
    }

    // Resize the pan/zoom drawings in `el`.
    // mode: undefined    -> resize only (keep the drawing exactly where it is,
    //                       so a tile drag only moves the controls);
    //       "keepcenter" -> keep the zoom but re-center the drawing within the
    //                       new bounds (window resize / fullscreen), so the view
    //                       stays consistent instead of drifting off-center.
    function resizePanZoomsIn(el, mode) {
        if (!el) return;
        el.querySelectorAll(".svg-host svg").forEach(function (svg) {
            if (svg._seismicPanZoom) resizePanZoom(svg._seismicPanZoom, mode);
        });
    }

    function resizePlotsIn(el, mode) {
        resizeFramedPlots(el);
        resizePanZoomsIn(el, mode);
    }

    // resize() moves the control icons back to the corner without touching the
    // viewport transform, so a plain resize leaves the drawing exactly where the
    // user put it (used for tile drags). "keepcenter" additionally re-centers the
    // drawing at the current zoom, using the container's new size, so the view is
    // consistently centered after a window resize or fullscreen (rather than
    // resetting the zoom or drifting off-center).
    function resizePanZoom(pz, mode) {
        try {
            pz.resize();
            if (mode === "keepcenter") { pz.center(); }
        } catch (e) { /* ignore */ }
    }
    function resizeAllPlots(mode) {
        Object.keys(grids).forEach(function (gid) {
            resizePlotsIn(grids[gid].el, mode);
        });
    }

    /* ---------------- SVG pan/zoom ---------------- */
    function initPanZoom() {
        if (typeof svgPanZoom === "undefined") return;
        document.querySelectorAll(".svg-host svg").forEach(function (svg) {
            try {
                var inst = svgPanZoom(svg, {
                    zoomEnabled: true,
                    controlIconsEnabled: true,
                    fit: true,
                    center: true,
                    minZoom: 0.3,
                    maxZoom: 20
                });
                svg._seismicPanZoom = inst;
                panZoomInstances.push(inst);
            } catch (e) { /* ignore */ }
        });
    }

    /* ---------------- GridStack per group ---------------- */
    var GRID_OPTS = {
        column: 12,
        cellHeight: 78,
        margin: 6,
        float: false,
        animate: true,
        handle: ".tile-header",
        resizable: { handles: "e, se, s, sw, w" },
        // scroll:false disables GridStack's own drag autoscroll, which keys off
        // the dragged tile's bounds (a tall plot scrolls constantly); we drive
        // autoscroll from the pointer position instead (see dragAutoScroll).
        draggable: { cancel: ".tile-btn", scroll: false }
    };

    // Pointer-position-based autoscroll while dragging a tile: the page scrolls
    // only when the cursor nears the top/bottom edge of the window, regardless
    // of how tall the dragged plot is.
    var _dragPointerY = null, _dragScrollRAF = null;
    function _onDragPointerMove(e) { _dragPointerY = e.clientY; }
    function _dragScrollTick() {
        if (_dragPointerY != null) {
            var margin = 64, maxSpeed = 22, h = window.innerHeight, y = _dragPointerY, dy = 0;
            if (y < margin) dy = -maxSpeed * (margin - y) / margin;
            else if (y > h - margin) dy = maxSpeed * (y - (h - margin)) / margin;
            if (dy) window.scrollBy(0, dy);
        }
        _dragScrollRAF = requestAnimationFrame(_dragScrollTick);
    }
    function startDragAutoScroll() {
        if (_dragScrollRAF) return;
        _dragPointerY = null;
        document.addEventListener("pointermove", _onDragPointerMove, true);
        _dragScrollRAF = requestAnimationFrame(_dragScrollTick);
    }
    function stopDragAutoScroll() {
        document.removeEventListener("pointermove", _onDragPointerMove, true);
        if (_dragScrollRAF) cancelAnimationFrame(_dragScrollRAF);
        _dragScrollRAF = null;
        _dragPointerY = null;
    }

    function setupGrid(gridEl, isReinit) {
        var gid = gridEl.getAttribute("data-group-id");
        var grid = GridStack.init(GRID_OPTS, gridEl);
        grids[gid] = grid;

        if (!isReinit) {
            // Capture the markup default layout before applying saved state.
            defaultLayouts[gid] = grid.save(false);
            var saved = load(LAYOUT_KEY(gid));
            if (saved && saved.length) {
                programmatic(function () {
                    try { grid.load(saved); } catch (e) { /* fall back to default */ }
                });
            }
        }

        var persist = debounce(function () { store(LAYOUT_KEY(gid), grid.save(false)); }, 250);
        grid.on("change added removed", function () { if (!suppressPersist) persist(); });
        grid.on("resizestart", function () { document.body.classList.add("gs-resizing"); });
        grid.on("dragstart", function () {
            document.body.classList.add("gs-dragging");
            startDragAutoScroll();
        });
        grid.on("resizestop dragstop", function (ev, el) {
            document.body.classList.remove("gs-resizing");
            document.body.classList.remove("gs-dragging");
            stopDragAutoScroll();
            resizePlotsIn(el);
        });
        grid.on("resize", function (ev, el) { resizePlotsIn(el); });

        gridEl.querySelectorAll(".grid-stack-item").forEach(observeTile);
        return grid;
    }

    function initGrids() {
        document.querySelectorAll(".grid-stack").forEach(function (gridEl) {
            setupGrid(gridEl, false);
        });
    }

    // Re-initialize a grid in place (keeping its DOM, so iframes are not
    // reloaded). Moving a group reparents its grid; if that happened while the
    // group was collapsed (zero height), GridStack loses the item pixel sizes,
    // so we rebuild it once the group is visible again.
    function reinitGrid(gridEl) {
        if (!gridEl || !gridEl.gridstack) return;
        try { gridEl.gridstack.destroy(false); } catch (e) { /* ignore */ }
        setupGrid(gridEl, true);
    }
    function markGridStale(gridEl) { if (gridEl) gridEl._needsReinit = true; }
    function reinitIfNeeded(gridEl) {
        // Only rebuild once the grid is actually visible (has real height).
        if (gridEl && gridEl._needsReinit && gridEl.offsetParent !== null) {
            reinitGrid(gridEl);
            gridEl._needsReinit = false;
        }
    }

    // Keep plots filling their tile as it changes size for any reason.
    function observeTile(item) {
        var body = item.querySelector(".tile-body");
        if (!body || typeof ResizeObserver === "undefined") return;
        if (item._seismicRO) { try { item._seismicRO.disconnect(); } catch (e) { /* ignore */ } }
        // Only the framed plots follow tile-body size here; pan/zoom drawings
        // are resized by the explicit handlers so their mode is correct.
        var ro = new ResizeObserver(debounce(function () { resizeFramedPlots(item); }, 60));
        item._seismicRO = ro;
        ro.observe(body);
    }

    /* ---------------- Tile controls ---------------- */
    function tileItem(btn) { return btn.closest(".grid-stack-item"); }

    function toggleCollapse(btn) {
        var item = tileItem(btn);
        var gid = item.closest(".grid-stack").getAttribute("data-group-id");
        var grid = grids[gid];
        var collapsed = item.classList.toggle("tile-collapsed");
        if (collapsed) {
            item.dataset.prevH = item.getAttribute("gs-h") || grid.getGridItems().length;
            var node = item.gridstackNode;
            item.dataset.prevH = node ? node.h : 4;
            grid.update(item, { h: 1, minH: 1, noResize: true });
        } else {
            var h = parseInt(item.dataset.prevH || "4", 10);
            grid.update(item, { h: h, minH: 1, noResize: false });
            resizePlotsIn(item);
        }
        setIcon(btn, collapsed ? "expand" : "collapse");
    }

    function toggleFullscreen(btn) {
        var item = tileItem(btn);
        activateTile(item);  // ensure its plot/iframe is rendered before enlarging
        var isFs = item.classList.toggle("gs-fullscreen");
        var backdrop = document.getElementById("fs-backdrop");
        if (isFs) {
            backdrop.hidden = false;
            document.body.style.overflow = "hidden";
        } else {
            backdrop.hidden = true;
            document.body.style.overflow = "";
        }
        setIcon(btn, isFs ? "shrink" : "expand-fs");
        // Keep the zoom but re-center within the new window/tile bounds, so the
        // framing is consistent across entering/leaving fullscreen. Plotly plots
        // fill regardless.
        requestAnimationFrame(function () { resizePlotsIn(item, "keepcenter"); });
        setTimeout(function () { resizePlotsIn(item, "keepcenter"); }, 80);
    }

    function exitFullscreen() {
        var item = document.querySelector(".grid-stack-item.gs-fullscreen");
        if (item) {
            item.classList.remove("gs-fullscreen");
            document.getElementById("fs-backdrop").hidden = true;
            document.body.style.overflow = "";
            var btn = item.querySelector('.tile-btn[data-act="fullscreen"]');
            if (btn) setIcon(btn, "expand-fs");
            requestAnimationFrame(function () { resizePlotsIn(item, "keepcenter"); });
            setTimeout(function () { resizePlotsIn(item, "keepcenter"); }, 80);
        }
    }

    /* ---------------- Group controls ---------------- */
    function anyGroupCollapsed() {
        return document.querySelector(".group.collapsed") != null;
    }
    function updateExpandAllLabel() {
        // The label always states what the next click will do: if anything is
        // collapsed the button expands all; otherwise it collapses all.
        var btn = document.getElementById("expand-all-btn");
        if (!btn) return;
        var lbl = btn.querySelector(".label");
        if (lbl) lbl.textContent = anyGroupCollapsed() ? "Expand all" : "Collapse all";
    }

    function setGroupCollapsed(section, collapsed, save, animate) {
        var body = section.querySelector(".group-body");
        var gridEl = section.querySelector(".grid-stack");
        if (collapsed) {
            if (animate === false) {
                section.classList.add("collapsed");
                body.style.height = "0px";
            } else {
                // Animate from the current height down to zero.
                body.style.transition = "none";
                body.style.height = body.scrollHeight + "px";
                body.offsetHeight; // reflow so the start height is applied
                body.style.transition = "";
                section.classList.add("collapsed");
                body.style.height = "0px";
            }
        } else {
            section.classList.remove("collapsed");
            // Make the grid visible and correctly sized (rebuild if a move left
            // it stale) before measuring the target height.
            body.style.transition = "none";
            body.style.height = "auto";
            reinitIfNeeded(gridEl);
            var target = body.scrollHeight;
            if (animate === false) {
                body.style.transition = "";
                body.style.height = "auto";
            } else {
                // Animate from zero up to the measured full height, then release
                // to auto so the group adapts to later content changes.
                body.style.height = "0px";
                body.offsetHeight; // reflow
                body.style.transition = "";
                body.style.height = target + "px";
                var onEnd = function (e) {
                    if (e.target !== body || e.propertyName !== "height") return;
                    body.style.height = "auto";
                    body.removeEventListener("transitionend", onEnd);
                };
                body.addEventListener("transitionend", onEnd);
            }
            requestAnimationFrame(function () { resizeAllPlots(); });
        }
        if (save) store(COLLAPSE_KEY(section.getAttribute("data-group-id")), collapsed);
        updateExpandAllLabel();
    }

    function initGroups() {
        document.querySelectorAll(".group").forEach(function (section) {
            var gid = section.getAttribute("data-group-id");
            if (load(COLLAPSE_KEY(gid))) setGroupCollapsed(section, true, false, false);
            var toggle = section.querySelector(".group-collapse");
            if (toggle) toggle.addEventListener("click", function () {
                setGroupCollapsed(section, !section.classList.contains("collapsed"), true);
            });
            var resetBtn = section.querySelector('[data-act="reset-group"]');
            if (resetBtn) resetBtn.addEventListener("click", function () {
                resetGroup(gid);
                toast("Layout reset");
            });
        });
    }

    function resetGroup(gid) {
        var grid = grids[gid];
        if (!grid) return;
        remove(LAYOUT_KEY(gid));
        // Un-collapse any collapsed tiles before restoring positions.
        grid.getGridItems().forEach(function (item) {
            if (item.classList.contains("tile-collapsed")) {
                item.classList.remove("tile-collapsed");
                grid.update(item, { noResize: false, minH: 1 });
                var btn = item.querySelector('.tile-btn[data-act="collapse"]');
                if (btn) setIcon(btn, "collapse");
            }
        });
        grid.load(defaultLayouts[gid]);
        resizeAllPlots();
    }

    function resetAll() {
        remove(ORDER_KEY);
        // Expand every group first so grids are visible when rebuilt.
        document.querySelectorAll(".group").forEach(function (section) {
            setGroupCollapsed(section, false, true);
        });
        // Restore the original group order. This reparents the sections, which
        // leaves their grids' items zero-height, so rebuild each grid afterward.
        var orderChanged = defaultGroupOrder &&
            currentGroupOrder().join(" ") !== defaultGroupOrder.join(" ");
        if (orderChanged) applyGroupOrder(defaultGroupOrder);
        Object.keys(grids).forEach(function (gid) {
            if (orderChanged) reinitGrid(grids[gid].el);
            resetGroup(gid);
        });
        requestAnimationFrame(function () { resizeAllPlots(); });
        toast("All layouts reset");
    }

    /* ---------------- Group ordering (drag to reorder) ---------------- */
    function currentGroupOrder() {
        return Array.prototype.map.call(
            document.querySelectorAll("main > .group"),
            function (s) { return s.getAttribute("data-group-id"); }
        );
    }
    function applyGroupOrder(order) {
        var main = document.querySelector("main");
        if (!main || !order) return;
        order.forEach(function (gid) {
            var sec = document.getElementById(gid);
            if (sec && sec.parentNode === main) main.appendChild(sec);
        });
    }
    function initGroupSort() {
        var main = document.querySelector("main");
        if (!main || typeof Sortable === "undefined") return;
        Sortable.create(main, {
            handle: ".group-drag-handle",
            draggable: ".group",
            animation: 150,
            ghostClass: "group-sortable-ghost",
            dragClass: "group-sortable-drag",
            onStart: function () { document.body.classList.add("gs-dragging"); },
            onEnd: function (evt) {
                document.body.classList.remove("gs-dragging");
                store(ORDER_KEY, currentGroupOrder());
                // Reordering reparents the moved section, which can leave its
                // grid's items zero-height (especially if it was collapsed).
                // Rebuild it now if visible, otherwise on its next expand.
                if (evt.item) {
                    var gridEl = evt.item.querySelector(".grid-stack");
                    markGridStale(gridEl);
                    reinitIfNeeded(gridEl);
                }
                requestAnimationFrame(function () { resizeAllPlots(); });
            }
        });
    }

    /* ---------------- Search filter ---------------- */
    function initSearch() {
        var input = document.getElementById("search-input");
        if (!input) return;
        input.addEventListener("input", debounce(function () {
            var q = input.value.trim().toLowerCase();
            document.querySelectorAll(".group").forEach(function (section) {
                var gid = section.getAttribute("data-group-id");
                var grid = grids[gid];
                var visible = 0;
                section.querySelectorAll(".grid-stack-item").forEach(function (item) {
                    var hay = (item.getAttribute("data-search") || "").toLowerCase();
                    var match = !q || hay.indexOf(q) !== -1;
                    item.hidden = !match;
                    if (match) visible++;
                });
                section.hidden = q && visible === 0;
                if (grid && visible) { try { grid.compact(); } catch (e) {} }
            });
            resizeAllPlots();
        }, 160));
    }

    /* ---------------- Expand / collapse all groups ---------------- */
    function initGlobalControls() {
        var resetBtn = document.getElementById("reset-all-btn");
        if (resetBtn) resetBtn.addEventListener("click", resetAll);


        var expandBtn = document.getElementById("expand-all-btn");
        if (expandBtn) expandBtn.addEventListener("click", function () {
            // Expand everything if anything is collapsed; otherwise collapse all.
            // (Matches the button label, which states the next action.)
            var expand = anyGroupCollapsed();
            document.querySelectorAll(".group").forEach(function (section) {
                setGroupCollapsed(section, !expand, true);
            });
        });

        // Tile action delegation.
        document.body.addEventListener("click", function (e) {
            var btn = e.target.closest(".tile-btn");
            if (!btn) return;
            var act = btn.getAttribute("data-act");
            if (act === "collapse") toggleCollapse(btn);
            else if (act === "fullscreen") toggleFullscreen(btn);
        });

        var backdrop = document.getElementById("fs-backdrop");
        if (backdrop) backdrop.addEventListener("click", exitFullscreen);
        document.addEventListener("keydown", function (e) {
            if (e.key === "Escape") exitFullscreen();
        });

        // Smooth-scroll group nav.
        document.querySelectorAll(".groupnav a").forEach(function (link) {
            link.addEventListener("click", function (e) {
                e.preventDefault();
                var t = document.getElementById(this.getAttribute("href").slice(1));
                if (t) t.scrollIntoView({ behavior: "smooth", block: "start" });
            });
        });
    }

    /* ---------------- Icons ---------------- */
    var ICONS = {
        collapse: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="18 15 12 9 6 15"/></svg>',
        expand: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="6 9 12 15 18 9"/></svg>',
        "expand-fs": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="15 3 21 3 21 9"/><polyline points="9 21 3 21 3 15"/><line x1="21" y1="3" x2="14" y2="10"/><line x1="3" y1="21" x2="10" y2="14"/></svg>',
        shrink: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="4 14 10 14 10 20"/><polyline points="20 10 14 10 14 4"/><line x1="14" y1="10" x2="21" y2="3"/><line x1="3" y1="21" x2="10" y2="14"/></svg>'
    };
    function setIcon(btn, name) { if (ICONS[name]) btn.innerHTML = ICONS[name]; }

    /* ---------------- Utilities ---------------- */
    function debounce(fn, ms) {
        var t;
        return function () {
            var ctx = this, args = arguments;
            clearTimeout(t);
            t = setTimeout(function () { fn.apply(ctx, args); }, ms);
        };
    }
    var toastTimer;
    function toast(msg) {
        var el = document.getElementById("toast");
        if (!el) return;
        el.textContent = msg;
        el.classList.add("show");
        clearTimeout(toastTimer);
        toastTimer = setTimeout(function () { el.classList.remove("show"); }, 1800);
    }

    /* ---------------- Auto-hiding top bar (header + group nav) ----------------
     * The top bar is fully visible at the top of the page. Once scrolled down
     * it slides out of the way, and reappears when the pointer moves to the
     * top edge of the window or hovers the bar itself. */
    function initTopbar() {
        var topbar = document.getElementById("topbar");
        if (!topbar) return;
        var peeking = false;
        var shown = true;
        var hideTimer = null;

        function setShown(s) {
            if (s === shown) return;
            shown = s;
            topbar.classList.toggle("hidden", !s);
            // Wrap so rAF's timestamp arg isn't read as `refit` (which would
            // re-fit/reset drawings just because the bar was revealed).
            if (s) requestAnimationFrame(function () { resizeAllPlots(); });
        }
        function update() {
            var atTop = window.scrollY < 8;
            setShown(atTop || peeking);
        }
        function peek(on) {
            if (on) {
                if (hideTimer) { clearTimeout(hideTimer); hideTimer = null; }
                peeking = true;
                update();
            } else {
                if (hideTimer) clearTimeout(hideTimer);
                hideTimer = setTimeout(function () { peeking = false; update(); }, 220);
            }
        }

        window.addEventListener("scroll", update, { passive: true });
        // A thin fixed strip at the top reveals the bar on hover; it sits above
        // graph iframes so it works even when the pointer is over a plot.
        var trigger = document.getElementById("topbar-trigger");
        if (trigger) {
            trigger.addEventListener("mouseenter", function () { peek(true); });
            trigger.addEventListener("mouseleave", function () { peek(false); });
        }
        // Secondary trigger for pointer motion not intercepted by an iframe.
        document.addEventListener("mousemove", function (e) {
            if (e.clientY <= 6) peek(true);
        }, { passive: true });
        // Keep it visible while hovering the bar; hide shortly after leaving.
        topbar.addEventListener("mouseenter", function () { peek(true); });
        topbar.addEventListener("mouseleave", function () { peek(false); });
        update();
    }

    /* ---------------- Boot ---------------- */
    function boot() {
        initTheme();
        // Apply a saved group order before grids/iframes are set up, so the
        // sections are reordered while iframes are still unloaded (avoids reloads).
        defaultGroupOrder = currentGroupOrder();
        applyGroupOrder(load(ORDER_KEY));
        initGrids();       // apply saved layout & size tiles
        initGroups();
        initPanZoom();
        initSearch();
        initGlobalControls();
        initGroupSort();
        initTopbar();
        updateExpandAllLabel();
        initLazyTiles();   // render/load tiles only as they near the viewport
        // Final settle: make sure every plot fills its tile.
        requestAnimationFrame(function () { resizeAllPlots(); });
        setTimeout(function () { resizeAllPlots(); }, 300);
        // On a window/page resize, re-fit drawings so they stay properly sized.
        window.addEventListener("resize", debounce(function () {
            resizeAllPlots("keepcenter");
        }, 150));
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", boot);
    } else {
        boot();
    }
})();
