/* window.vala
 *
 * Copyright 2020 Ahmed Eldemery
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE X CONSORTIUM BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Except as contained in this notice, the name(s) of the above copyright
 * holders shall not be used in advertising or otherwise to promote the sale,
 * use or other dealings in this Software without prior written
 * authorization.
 */

namespace Gtk4LayoutmanagerTransformation {
    public class MainWindow : Gtk.ApplicationWindow {
        const string[] icon_names = {
            "action-unavailable-symbolic",
            "address-book-new-symbolic",
            "application-exit-symbolic",
            "appointment-new-symbolic",
            "bookmark-new-symbolic",
            "call-start-symbolic",
            "call-stop-symbolic",
            "camera-switch-symbolic",
            "chat-message-new-symbolic",
            "color-select-symbolic",
            "contact-new-symbolic",
            "document-edit-symbolic",
            "document-new-symbolic",
            "document-open-recent-symbolic",
            "document-open-symbolic",
            "document-page-setup-symbolic",
            "document-print-preview-symbolic",
            "document-print-symbolic",
            "document-properties-symbolic",
            "document-revert-symbolic-rtl",
            "document-revert-symbolic",
            "document-save-as-symbolic",
            "document-save-symbolic",
            "document-send-symbolic",
            "edit-clear-all-symbolic",
            "edit-clear-symbolic-rtl",
            "edit-clear-symbolic",
            "edit-copy-symbolic",
            "edit-cut-symbolic",
            "edit-delete-symbolic",
            "edit-find-replace-symbolic",
            "edit-find-symbolic",
            "edit-paste-symbolic",
            "edit-redo-symbolic-rtl",
            "edit-redo-symbolic",
            "edit-select-all-symbolic",
            "edit-select-symbolic",
            "edit-undo-symbolic-rtl",
            "edit-undo-symbolic",
            "error-correct-symbolic",
            "find-location-symbolic",
            "folder-new-symbolic",
            "font-select-symbolic",
            "format-indent-less-symbolic-rtl",
            "format-indent-less-symbolic",
            "format-indent-more-symbolic-rtl",
            "format-indent-more-symbolic",
            "format-justify-center-symbolic",
            "format-justify-fill-symbolic",
            "format-justify-left-symbolic",
            "format-justify-right-symbolic",
            "format-text-bold-symbolic",
            "format-text-direction-symbolic-rtl",
            "format-text-direction-symbolic",
            "format-text-italic-symbolic",
            "format-text-strikethrough-symbolic",
            "format-text-underline-symbolic",
            "go-bottom-symbolic",
            "go-down-symbolic",
            "go-first-symbolic-rtl",
            "go-first-symbolic",
            "go-home-symbolic",
            "go-jump-symbolic-rtl",
            "go-jump-symbolic",
            "go-last-symbolic-rtl",
            "go-last-symbolic",
            "go-next-symbolic-rtl",
            "go-next-symbolic",
            "go-previous-symbolic-rtl",
            "go-previous-symbolic",
            "go-top-symbolic",
            "go-up-symbolic",
            "help-about-symbolic",
            "insert-image-symbolic",
            "insert-link-symbolic",
            "insert-object-symbolic",
            "insert-text-symbolic",
            "list-add-symbolic",
            "list-remove-all-symbolic",
            "list-remove-symbolic",
            "mail-forward-symbolic",
            "mail-mark-important-symbolic",
            "mail-mark-junk-symbolic",
            "mail-mark-notjunk-symbolic",
            "mail-message-new-symbolic",
            "mail-reply-all-symbolic",
            "mail-reply-sender-symbolic",
            "mail-send-receive-symbolic",
            "mail-send-symbolic",
            "mark-location-symbolic",
            "media-eject-symbolic",
            "media-playback-pause-symbolic",
            "media-playback-start-symbolic",
            "media-playback-stop-symbolic",
            "media-record-symbolic",
            "media-seek-backward-symbolic",
            "media-seek-forward-symbolic",
            "media-skip-backward-symbolic",
            "media-skip-forward-symbolic",
            "media-view-subtitles-symbolic",
            "object-flip-horizontal-symbolic",
            "object-flip-vertical-symbolic",
            "object-rotate-left-symbolic",
            "object-rotate-right-symbolic",
            "object-select-symbolic",
            "open-menu-symbolic",
            "process-stop-symbolic",
            "send-to-symbolic",
            "sidebar-hide-symbolic",
            "sidebar-show-symbolic",
            "star-new-symbolic",
            "system-log-out-symbolic",
            "system-reboot-symbolic",
            "system-run-symbolic",
            "system-search-symbolic",
            "system-shutdown-symbolic",
            "system-switch-user-symbolic",
            "tab-new-symbolic",
            "tools-check-spelling-symbolic",
            "value-decrease-symbolic",
            "value-increase-symbolic",
            "view-app-grid-symbolic",
            "view-conceal-symbolic",
            "view-continuous-symbolic",
            "view-dual-symbolic",
            "view-fullscreen-symbolic",
            "view-grid-symbolic",
            "view-list-bullet-symbolic",
            "view-list-ordered-symbolic",
            "view-list-symbolic",
            "view-mirror-symbolic",
            "view-more-horizontal-symbolic",
            "view-more-symbolic",
            "view-paged-symbolic",
            "view-pin-symbolic",
            "view-refresh-symbolic",
            "view-restore-symbolic",
            "view-reveal-symbolic",
            "view-sort-ascending-symbolic",
            "view-sort-descending-symbolic",
            "zoom-fit-best-symbolic",
            "zoom-in-symbolic",
            "zoom-original-symbolic",
            "zoom-out-symbolic",
        };
        
        public MainWindow (Gtk.Application app) {
            Object (application: app);
            title = "Layout Manager - Transformation";
            default_width = 600;
            default_height = 620;

            var globe_widget = new GlobeWidget ();

            foreach (var name in icon_names) {
                var child = new Gtk.Image.from_icon_name (name);
                child.margin_start = child.margin_end = child.margin_top = child.margin_bottom = 4;

                globe_widget.add_child (child);
            }
            
            this.set_child (globe_widget);
        }
    }
}
