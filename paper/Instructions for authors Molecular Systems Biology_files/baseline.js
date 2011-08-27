(function($) {
	
	$.fn.baseline = function(childSelector) {
		if (!childSelector) {
			childSelector = '.baseline';
		}
		var self = this;
		var getElementInfo = function(el) {
			var offset = $(el).offset();
			var width = $(el).outerWidth();
			var height = $(el).outerHeight();
			return {
				el: el,
				width: width,
				height: height,
				left: offset.left,
				right: offset.left + width,
				top: offset.top,
				bottom: offset.top + height
			};
		};
		var groupItems = function(items, mode) {
			var m = -1, n = items.length;
			if (!n) {
				return false;
			}
			
			if (mode == 'rows') {
				var testProp = 'top';
				var groupProp = 'bottom';
			} else {
				var testProp = 'left';
				var groupProp = 'right';
			}
			
			var group = [];
			for (var i = 0; i < n; ++i) {
				if (m == -1 || items[i][testProp] > group[m].max) {
					group[++m] = {max: items[i][groupProp], items: [items[i]]};
				} else {
					group[m].max = Math.max(group[m].max, items[i][groupProp]);
					group[m].items.push(items[i]);
				}
			}
			return group;
		};
		
		var setAlignment = function() {
			var blocks = (function() {
				// for each of the wrapper elements group together the children for resizing
				var b = [], i = -1;
				self.each(function() {
					b[++i] = [];
					$(this).find(childSelector).each(function() {
						b[i].push(getElementInfo(this));
					});
				});
				return b;
			})();
			if (!blocks.length) {
				return;
			}
			
			$.each(blocks, function(index, items) {
				items.sort(function(a, b) {
					return a.bottom > b.bottom ? 1 : -1;
				});
				
				var rows = groupItems(items, 'rows');
				for (var i = 0; i < rows.length; ++i) {
					rowItems = rows[i].items;
					rowItems.sort(function(a, b) {
						return a.right > b.right ? 1 : -1;
					});
					// within each row, group items into columns
					var cols = groupItems(rowItems, 'cols');
					for (var j = 0; j < cols.length; ++j) {
						var colItems = cols[j].items;
						colItems.sort(function(a, b) {
							return a.bottom > b.bottom ? 1 : -1;
						});
						// if there is more than 1 item in a sub column, only adjust the bottom	item
						var item = colItems.pop();
						var d = rows[i].max - item.bottom;
						if (d > 0) {
							var h = !$.boxModel ? $(item.el).outerHeight() : $(item.el).height();
							$(item.el).height(h + d);
						}
					}
				}
			});
		};
		$(window).bind('resize', setAlignment);
		setTimeout(function() {
			setAlignment();
		}, 100);  // timeout to deal with IE JS minmax resizing
		return this;	
	};
})(jQuery);
