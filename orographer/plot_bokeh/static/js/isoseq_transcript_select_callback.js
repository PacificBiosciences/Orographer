export default function (args) {
const transcript_select = args.transcript_select;
const transcript_sources = args.transcript_sources;
const selectedTranscript = transcript_select.value || "ALL";

function setSelectedIndices(source, indices) {
    source.selected.indices = indices;
    source.selected.change.emit();
}

transcript_sources.forEach(function (source) {
    const data = source.data;
    const matches = [];
    if (selectedTranscript != "ALL") {
        (data.transcript_id || []).forEach(function (transcriptId, rowIndex) {
            if (String(transcriptId) == String(selectedTranscript)) {
                matches.push(rowIndex);
            }
        });
    }
    setSelectedIndices(source, matches);
    source.change.emit();
});
}
