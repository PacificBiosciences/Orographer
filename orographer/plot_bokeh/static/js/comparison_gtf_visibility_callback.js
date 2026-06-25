export default function (args) {
const hideComparisonCheckbox = args.hide_comparison_checkbox;
const comparisonComponents = args.comparison_components || [];
const isVisible = hideComparisonCheckbox.active ? false : true;

comparisonComponents.forEach(function (component) {
    if (!component) return;
    component.visible = isVisible;
    component.change.emit();
});
}
