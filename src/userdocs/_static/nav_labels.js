/*
 * Insert group-heading labels into the "How to Use" sidebar section.
 *
 * Sphinx/Furo strips toctree :caption: labels from the sidebar navigation
 * tree, so we re-add them by injecting <div class="nav-group-label"> elements
 * before the first item of each group.  Using a.pathname (the browser-resolved
 * absolute path) means the match works regardless of whether the href is a
 * relative path or "#" (current page).
 */
document.addEventListener("DOMContentLoaded", function () {
    var sidebar = document.querySelector(".sidebar-tree");
    if (!sidebar) return;

    // Each entry: path suffix that identifies the first item in the group,
    // and the label text to inject above it.
    var groups = [
        { suffix: "/workflow/wf.html",        label: "Workflow commands" },
        { suffix: "/utility/splitbam.html",   label: "Utility commands"  },
        { suffix: "/use/inputs.html",         label: "Shared topics"     },
    ];

    sidebar.querySelectorAll("li.toctree-l2 > a").forEach(function (a) {
        var pathname = a.pathname || "";
        groups.forEach(function (group) {
            if (pathname.endsWith(group.suffix)) {
                var div = document.createElement("div");
                div.className = "nav-group-label";
                div.textContent = group.label;
                a.closest("li.toctree-l2").before(div);
            }
        });
    });
});
